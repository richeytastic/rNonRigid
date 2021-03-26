/************************************************************************
 * Copyright (C) 2021 Richard Palmer
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
#include <cassert>
using rNonRigid::ViscoElasticTransformer;
using rNonRigid::SmoothingWeights;
using rNonRigid::MatX6f;
using rNonRigid::MatX3f;
using rNonRigid::MatX3f;
using rNonRigid::MatXf;
using rNonRigid::MatXi;
using rNonRigid::VecXf;
using rNonRigid::Vec3f;


namespace {

MatXf createWeights( const SmoothingWeights &swts, const VecXf &iwts)
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


void regularise( MatX3f &M, const MatXf &wts, const VecXf &wRowSums, const MatXi &nidxs, size_t nSteps)
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


// Outlier diffusion looks at relatively low probability vectors and averages them with their
// neighbours to diffuse possible error among the neighbours. Higher probabilities are ignored.
void diffuseOutliers( MatX3f &M,
                      const SmoothingWeights &swts, const VecXf &iwts, float wthresh, size_t nSteps)
{
    // Identify outliers as those with weights lower than threshold
    int N = int(M.rows());
    std::vector<int> outliers;
    outliers.reserve(N);
    for ( int i = 0; i < N; ++i)
        if ( iwts[i] < wthresh)
            outliers.push_back(i);
    N = int(outliers.size());

    const size_t K = swts.indices().cols();

    // Calculate sum of weights over all neighbours of each outlier
    VecXf wsums = VecXf::Zero( N);
    for ( int i = 0; i < N; ++i)
    {
        const int l = outliers[i];
        for ( size_t k = 0; k < K; ++k)
            wsums[i] += swts.weights()(l,k);
    }   // end for

    for ( size_t s = 0; s < nSteps; ++s)
    {
        // Making a copy of all vertices each time so that neighbour effects
        // aren't compounded (low threshold vertices having low threshold neighbours
        // and vice versa) does not empirically result in an improved registration
        // so this is removed.
        //const MatX3f tmpM = M;

        for ( int i = 0; i < N; ++i)
        {
            const int l = outliers[i];
            Vec3f vavg = Vec3f::Zero();
            for ( size_t k = 0; k < K; ++k)
            {
                const float wt = swts.weights()(l,k);
                //vavg += wt * tmpM.row( swts.indices()(l,k));
                vavg += wt * M.row( swts.indices()(l,k));
            }   // end for

            // The amount of the neighbouring vector directions added is in
            // proportion to the uncertainty about the direction of this vector.
            const float iw = iwts[l];
            M.row(l) *= iw;
            M.row(l) += (1.0f - iw) * vavg / wsums[i];
        }   // end for
    }   // end for
}   // end diffuseOutliers

}  // end namespace


ViscoElasticTransformer::ViscoElasticTransformer( const SmoothingWeights &swts,
                                                  size_t nvs, size_t nve,
                                                  size_t nes, size_t nee,
                                                  size_t numUpdates, float itw, size_t nodi)
    : _swts(swts),
      _viscousAnnealingRate( std::exp( std::log( float(nve)/float(nvs)) / numUpdates)),
      _elasticAnnealingRate( std::exp( std::log( float(nee)/float(nes)) / numUpdates)),
      _numViscousStart( nvs),
      _numElasticStart( nes),
      _inlierThresholdWt( itw),
      _numOutlierDiffIts(nodi),
      _field( MatX3f::Zero( swts.indices().rows(), 3)), // total displacement field
      _i(0.0f)
{
}   // end ctor


void ViscoElasticTransformer::update( MatX3f &df, const VecXf &iwts)
{
    assert( df.rows() == int(_swts.indices().rows()));
    assert( df.rows() == iwts.size());

    const MatXf wts = createWeights( _swts, iwts);
    const VecXf wRowSums = wts.rowwise().sum();

    // Regularise given displacement field prior to adding to total displacement
    const size_t nVs = size_t( _numViscousStart * std::pow( _viscousAnnealingRate, _i));
    regularise( df, wts, wRowSums, _swts.indices(), nVs);
    const MatX3f pfield = _field;   // Copy prior field
    _field += df;

    // Regularise/relax the TOTAL deformation field
    const size_t nEs = size_t( _numElasticStart * std::pow( _elasticAnnealingRate, _i));
    regularise( _field, wts, wRowSums, _swts.indices(), nEs);
    diffuseOutliers( _field, _swts, iwts, _inlierThresholdWt, _numOutlierDiffIts);
    _i += 1.0f;

    df = _field - pfield;   // Set the difference in the deformation field
}   // end update

