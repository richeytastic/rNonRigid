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

#include <InlierFinder.h>
#include <cassert>
#include <cfloat>
#include <cmath>
using rNonRigid::InlierFinder;
using rNonRigid::MatXf;
using rNonRigid::VecXf;


InlierFinder::InlierFinder( float k, bool uo, size_t n)
    : _kappa(k), _useOrientation(uo), _numIterations(n), _minSig( 1.0f/n), _maxSig( (float)n) {}


VecXf InlierFinder::operator()( const MatXf &rfA, const MatXf &rfB, const VecXf &flags) const
{
    const size_t N = rfA.rows();
    assert( long(N) == rfB.rows());
    assert( long(N) == flags.size());

    VecXf probs = flags;

    // Calculate distance squared deltas over the features. In the original implementation,
    // these values are repeatedly calculated within the _numIterations loop which is inefficient
    // since the values are not being updated.
    //
    // In the original implementation the squared distances are calculated over the whole (6D) feature.
    // However, I am unsure if this was actually the intention since the orientation component is
    // explicitly dealt with separately after the main iteration loop. After having tested taking
    // the difference over the whole feature versus just the position component, the empirical
    // difference is virtually impossible to detect so the more efficient computation is used here.
    const VecXf l2sqs = (rfB.leftCols<3>() - rfA.leftCols<3>()).rowwise().squaredNorm();

    static const float G_CONST = 1.0f/std::sqrt(float(2.0 * EIGEN_PI));
    const float G_FACTOR = G_CONST * std::exp( -0.5f * _kappa * _kappa);

    // The number of iterations is given as a magic number of 10 in the MeshMonk implementation.
    // This value is not justified in the implementation or the supplementary material so is
    // included here as a parameter for testing.
    for ( size_t i = 0; i < _numIterations; ++i)
    {
        const float sigDen = probs.sum() + FLT_MIN;
        assert( !std::isnan(sigDen));
        const float sigNum = probs.dot(l2sqs);
        const float ssnd = std::sqrt(sigNum/sigDen);
        assert( !std::isnan(ssnd));

        // Recalc position based probabilities
        const float sigma = std::max( _minSig, std::min( ssnd, _maxSig));
        const float expFact = -0.5f / (sigma * sigma);
        assert( !std::isnan(expFact));
        const float lambda = G_FACTOR / sigma;
        const float gConst = G_CONST / sigma;

        for ( size_t j = 0; j < N; ++j)
        {
            const float p = gConst * std::exp( expFact * l2sqs[j]);
            probs[j] *= p / (p + lambda);
        }   // end for
    }   // end for

    // Decrease weights the less congruent the surface normals are?
    if ( _useOrientation)
    {
        static const float EPS = 1e-6f;
        static const float ONE_MINUS_EPS = 1.0f - EPS;
        // Scale the dot products of the respective normals to be in [EPS, 1.0f]
        const VecXf d = ONE_MINUS_EPS * (0.5f * (rfA.rightCols<3>().array() * rfB.rightCols<3>().array()).rowwise().sum() + 0.5f) + EPS;
        //if ( (d.sum() / N) < 0.5f)
        //    std::cerr << "[WARNING] rNonRigid::InlierFinder: Very low inlier weights due to surface normals." << std::endl;
        probs = probs.array() * d.array();
    }   // end if

    return probs;
}   // end operator()
