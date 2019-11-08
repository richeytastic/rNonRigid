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

#include <InlierFinder.h>
#include <iostream>
#include <cassert>
#include <cmath>
using RNonRigid::InlierFinder;
using RNonRigid::VecXf;
using RNonRigid::FlagVec;
using RNonRigid::FeatMat;


InlierFinder::InlierFinder( float k, bool uo, size_t n)
    : _kappa(k), _useOrientation(uo), _numIterations(n), _minSig( 1.0f/n), _maxSig( (float)n) {}


VecXf InlierFinder::operator()( const FeatMat &rfA, const FeatMat &rfB, const FlagVec &flags) const
{
    const size_t N = rfA.rows();
    assert( rfA.rows() == rfB.rows());
    assert( rfA.rows() == flags.size());
    assert( rfA.cols() == rfB.cols());

    VecXf probs = flags;

    // Calculate distance squared deltas over the 6D features. In the original implementation,
    // these values are repeatedly calculated within the _numIterations loop which is inefficient
    // since the values are not being updated.
    //
    // In the original implementation the squared distances are calculated over the whole 6D feature.
    // This is replicated here. However, I am unsure if this was actually the intention since the
    // orientation component is explicitly dealt with separately after the main iteration loop.
    // Probably doesn't make a huge practical difference, but needs justifying nevertheless.
    std::vector<float> l2sqs(N);
    for ( size_t i = 0; i < N; ++i)
        l2sqs[i] = (rfB.row(i) - rfA.row(i)).squaredNorm();

    static const double G_CONST = 1.0/sqrt(2.0 * EIGEN_PI);
    const double G_FACTOR = G_CONST * exp( -0.5 * _kappa * _kappa);

    // The number of iterations is given as a magic number of 10 in the MeshMonk implementation.
    // This value is not justified in the implementation or the supplementary material so is
    // included here as a parameter for testing.
    for ( size_t i = 0; i < _numIterations; ++i)
    {
        const double sigDen = probs.sum();
        assert( sigDen > 0.0);
        double sigNum = 0.0;
        for ( size_t j = 0; j < N; ++j)
            sigNum += probs[j] * l2sqs[j];

        const double sigma = std::max<double>( _minSig, std::min<double>( sqrt( sigNum/sigDen), _maxSig));
        const double lambda = G_FACTOR / sigma;

        // Recalc distance based probabilities
        const double sigmaSq = sigma * sigma;
        for ( size_t j = 0; j < N; ++j)
        {
            const double p = G_CONST * exp( -0.5 * l2sqs[j] / sigmaSq) / sigma;
            probs[j] *= float(p / (p + lambda));
        }   // end for
    }   // end for

    if ( _useOrientation)
    {
        float avgOrientationInlierWeight = 0.0f;    // For warning when too low (may mean normals are flipped)
        for ( size_t i = 0; i < N; ++i)
        {
            const float dp = rfA.row(i).tail<3>().dot( rfB.row(i).tail<3>());
            const float p = std::min( 1.0f, 0.5f + dp / 2.0f); // Rescale to be in [0,1]
            avgOrientationInlierWeight += p;
            probs[i] *= p;
        }   // end for

        avgOrientationInlierWeight /= N;
        if ( avgOrientationInlierWeight < 0.5f)
            std::cerr << "[WARNING] RNonRigid::InlierFinder: Very low inlier weights due to surface normals." << std::endl;
    }   // end if

    return probs;
}   // end operator()
