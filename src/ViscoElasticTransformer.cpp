/************************************************************************
 * Copyright (C) 2019 Richard Palmer
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
#include <cassert>
using RNonRigid::ViscoElasticTransformer;
using RNonRigid::SmoothingWeights;
using RNonRigid::FeatMat;
using RNonRigid::DispMat;
using RNonRigid::FlagVec;
using RNonRigid::MatXf;
using RNonRigid::MatXi;
using RNonRigid::VecXf;
using RNonRigid::Vec3f;


namespace {

// Regularise the given displacement map
void regulariseField( DispMat &DF, const SmoothingWeights &Swts, const VecXf &iws, size_t nSteps)
{
    const float EPS = 1e-5f;
    const MatXi& nidxs = Swts.indices();    // Neighbour indices
    const MatXf& swts = Swts.weights();     // Smoothing weights

    const size_t K = nidxs.cols();  // Number of neighbours of each vertex to iterate over
    const size_t N = DF.rows();     // Number of vertices in field

    DispMat rDF( N, 3);

    for ( size_t it = 0; it < nSteps; ++it)
    {
        for ( size_t i = 0; i < N; ++i)
        {
            Vec3f vavg = Vec3f::Zero();
            float wsum = 0.0f;

            for ( size_t k = 0; k < K; ++k)
            {
                const size_t j = nidxs(i,k);   // Neighbour index
                const float wt = (1.0f - EPS) * (iws[j] * swts(i,k)) + EPS;  // Rescale weight in [EPS,1]
                vavg += DF.row(j) * wt;
                wsum += wt;
            }   // end for

            rDF.row(i) = vavg / wsum;
        }   // end for

        DF = rDF;
    }   // end for
}   // end regulariseField


// Outlier diffusion looks at relatively low probability vectors and averages these with their neighbours
// to diffuse the possible error among its neighbours. Higher probability vectors are left alone.
void diffuseOutliers( DispMat &DF, const SmoothingWeights &Swts, const VecXf &iws, float wthresh, size_t nSteps)
{
    const MatXi& nidxs = Swts.indices();    // Neighbour indices
    const MatXf& swts = Swts.weights();     // Smoothing weights

    // Identify outliers outside of iterative loop for efficiency
    int N = int(DF.rows());
    std::vector<int> outliers;
    for ( int i = 0; i < N; ++i)
    {
        if ( iws[i] <= wthresh)
            outliers.push_back(i);
    }   // end for
    N = int(outliers.size());

    const size_t K = nidxs.cols();
    for ( size_t s = 0; s < nSteps; ++s)
    {
        const DispMat tmpDF = DF;

        for ( int i = 0; i < N; ++i)
        {
            const int l = outliers[i];
            const float iw = iws[l];

            float wsum = 0.0f;
            Vec3f vavg = Vec3f::Zero();
            for ( size_t k = 0; k < K; ++k)
            {
                const float wt = swts(l,k);
                vavg += wt * tmpDF.row( nidxs(l,k));
                wsum += wt;
            }   // end for

            // The amount of the neighbouring vector directions added is in
            // proportion to the uncertainty about the direction of this vector.
            DF.row(l) *= iw;
            DF.row(l) += (1.0f - iw) * vavg / wsum;
        }   // end for
    }   // end for
}   // end diffuseOutliers


}   // end namespace



ViscoElasticTransformer::ViscoElasticTransformer( size_t nvs, size_t nve, size_t nes, size_t nee, size_t numUpdates, float itw, size_t nodi)
    : _viscousAnnealingRate( expf( logf( float(nve)/float(nvs)) / (numUpdates-1))),
      _elasticAnnealingRate( expf( logf( float(nee)/float(nes)) / (numUpdates-1))),
      _numViscousStart( nvs),
      _numElasticStart( nes),
      _inlierThresholdWt( itw),
      _numOutlierDiffIts(nodi)
{
#ifndef NDEBUG
    std::cout << "Viscous rate: " << _viscousAnnealingRate << " | start: " << nvs << " | end: " << nve << std::endl;
    std::cout << "Elastic rate: " << _elasticAnnealingRate << " | start: " << nes << " | end: " << nee << std::endl;
    std::cout << "Iterations: " << numUpdates << std::endl;
#endif 
}   // end ctor



DispMat ViscoElasticTransformer::update( size_t i, const DispMat &CF, const SmoothingWeights &swts, const VecXf &iws, const FlagVec &flgs)
{
    assert( CF.rows() == int(swts.indices().rows()));
    assert( CF.rows() == iws.size());
    assert( CF.rows() == flgs.size());

    if ( i == 0)    // Reset the displacement field?
        _thisDispField = DispMat::Zero( CF.rows(), 3);

    const size_t numViscousIts = size_t( roundf( _numViscousStart * powf( _viscousAnnealingRate, float(i))));
    const size_t numElasticIts = size_t( roundf( _numElasticStart * powf( _elasticAnnealingRate, float(i))));

    _prevDispField = _thisDispField;                  // Remember from previous iteration
    DispMat rDF = CF.array().colwise() * iws.array(); // Component wise weighting by the inliers
    regulariseField( rDF, swts, iws, numViscousIts);  // Regularise the displacement field
    _thisDispField += rDF;                            // Viscous addition

    regulariseField( _thisDispField, swts, iws, numElasticIts); // Regularise TOTAL deformation field for elastic effect

    diffuseOutliers( _thisDispField, swts, iws, _inlierThresholdWt, _numOutlierDiffIts);

    return _thisDispField - _prevDispField;
}   // end update


