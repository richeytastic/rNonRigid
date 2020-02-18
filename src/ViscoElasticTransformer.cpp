/************************************************************************
 * Copyright (C) 2020 Richard Palmer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/

#include <ViscoElasticTransformer.h>
#include <iostream>
#include <iomanip>
#include <cassert>
using rNonRigid::ViscoElasticTransformer;
using rNonRigid::SmoothingWeights;
using rNonRigid::FeatMat;
using rNonRigid::MatX3f;
using rNonRigid::MatX3f;
using rNonRigid::FlagVec;
using rNonRigid::MatXf;
using rNonRigid::MatXi;
using rNonRigid::VecXf;
using rNonRigid::Vec3f;


MatXf rNonRigid::createWeights( const SmoothingWeights &swts, const VecXf &iwts)
{
    static const float EPS = 1e-5f;
    static const float ONE_MINUS_EPS = 1.0f - EPS;
    const size_t N = swts.indices().rows();  // Number of vertices in field
    const size_t K = swts.indices().cols();  // Number of neighbours of each vertex to iterate over

    MatXf wts(N, K);
    for ( size_t i = 0; i < N; ++i)
    {
        for ( size_t k = 0; k < K; ++k) // Typically 80 or so vertices nearest to i
        {
            const size_t j = swts.indices()(i,k);   // Neighbour index
            wts(i,k) = ONE_MINUS_EPS * iwts[j] * swts.weights()(i,k) + EPS;  // Rescale weight in [EPS,1]
        }   // end for
    }   // end for

    return wts;
}   // end createWeights


void rNonRigid::regularise( MatX3f &M, const MatXf &wts, const VecXf &wRowSums, const MatXi &nidxs, size_t nSteps)
{
    const size_t N = wts.rows();  // Number of vertices in field
    const size_t K = wts.cols() - 1;  // Number of neighbours of each vertex to iterate over
    MatX3f rM( N, 3);

    for ( size_t it = 0; it < nSteps; ++it)
    {
        for ( size_t i = 0; i < N; ++i)
        {
            Vec3f vavg = Vec3f::Zero();
            for ( size_t k = 0; k < K; ++k) // Typically 80 or so vertices nearest to i
                vavg += M.row(nidxs(i,k)) * wts(i,k);
            rM.row(i) = vavg / wRowSums[i];
        }   // end for

        M = rM;
    }   // end for
}   // end regularise


// Outlier diffusion looks at relatively low probability vectors and averages these with their neighbours
// to diffuse the possible error among its neighbours. Higher probability vectors are left alone.
void rNonRigid::diffuseOutliers( MatX3f &M, const SmoothingWeights &swts, const VecXf &iwts, float wthresh, size_t nSteps)
{
    // Identify outliers outside of iterative loop for efficiency
    int N = int(M.rows());
    std::vector<int> outliers;
    for ( int i = 0; i < N; ++i)
    {
        if ( iwts[i] <= wthresh)
            outliers.push_back(i);
    }   // end for
    N = int(outliers.size());

    const size_t K = swts.indices().cols();
    for ( size_t s = 0; s < nSteps; ++s)
    {
        const MatX3f tmpM = M;

        for ( int i = 0; i < N; ++i)
        {
            const int l = outliers[i];
            const float iw = iwts[l];

            float wsum = 0.0f;
            Vec3f vavg = Vec3f::Zero();
            for ( size_t k = 0; k < K; ++k)
            {
                const float wt = swts.weights()(l,k);
                vavg += wt * tmpM.row( swts.indices()(l,k));
                wsum += wt;
            }   // end for

            // The amount of the neighbouring vector directions added is in
            // proportion to the uncertainty about the direction of this vector.
            M.row(l) *= iw;
            M.row(l) += (1.0f - iw) * vavg / wsum;
        }   // end for
    }   // end for
}   // end diffuseOutliers


ViscoElasticTransformer::ViscoElasticTransformer( size_t nvs, size_t nve, size_t nes, size_t nee, size_t numUpdates, float itw, size_t nodi)
    : _viscousAnnealingRate( expf( logf( float(nve)/float(nvs)) / (numUpdates-1))),
      _elasticAnnealingRate( expf( logf( float(nee)/float(nes)) / (numUpdates-1))),
      _numViscousStart( nvs),
      _numElasticStart( nes),
      _inlierThresholdWt( itw),
      _numOutlierDiffIts(nodi)
{
    //std::cout << "ElasticAnnealingRate = " << _elasticAnnealingRate << std::endl;
}   // end ctor


void ViscoElasticTransformer::update( size_t i, MatX3f &D, const SmoothingWeights &swts, const VecXf &iwts)
{
    assert( D.rows() == int(swts.indices().rows()));
    assert( D.rows() == iwts.size());

    if ( _thisField.size() == 0)    // Reset the displacement field?
        _thisField = MatX3f::Zero( D.rows(), 3);

    const MatXf wts = createWeights( swts, iwts);
    const VecXf wRowSums = wts.rowwise().sum();

    const size_t numViscousIts = size_t( roundf( _numViscousStart * powf( _viscousAnnealingRate, float(i))));
    const size_t numElasticIts = size_t( roundf( _numElasticStart * powf( _elasticAnnealingRate, float(i))));
    //std::cout << std::setw(3) << i << ") numElasticIts = " << numElasticIts << std::endl;

    _prevField = _thisField;                                        // Remember from previous iteration
    regularise( D, wts, wRowSums, swts.indices(), numViscousIts);   // Regularise the displacement field
    _thisField += D;                                                // Viscous addition

    // Regularise TOTAL deformation field for elastic effect
    regularise( _thisField, wts, wRowSums, swts.indices(), numElasticIts);

    diffuseOutliers( _thisField, swts, iwts, _inlierThresholdWt, _numOutlierDiffIts);

    D = _thisField - _prevField;
}   // end update


