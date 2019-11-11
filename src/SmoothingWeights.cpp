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

#include <SmoothingWeights.h>
#include <cassert>
#include <cmath>
using rNonRigid::SmoothingWeights;
using rNonRigid::KNNMap;
using rNonRigid::FlagVec;


SmoothingWeights::SmoothingWeights( const KNNMap &kmap, const FlagVec &flgs, float sigma)
    : _kmap(kmap), _smw( kmap.indices().rows(), kmap.indices().cols())
{
    static const float EPS = 1e-5f;
    const float EXP_FACTOR = -0.5f / powf( sigma, 2.0f);
    const size_t K = _smw.cols();
    const size_t N = _smw.rows();

    for ( size_t i = 0; i < N; ++i)
    {
        float wsum = 0.0f;
        for ( size_t k = 0; k < K; ++k)
        {
            const float dsq = kmap.sqDiffs()(i,k);      // Squared distance
            const float gwt = expf( dsq * EXP_FACTOR);  // Gaussian weight
            const int j = kmap.indices()(i,k);          // Neighbour
            const float nflg = flgs[j];                 // Neighbour flag

            // Rescale weight in [EPS,1] instead of [0,1] so that nodes with an inlier weight
            // of 0 don't end up with a deformation vector having zero magnitude.
            const float wt = (1.0f - EPS) * (nflg * gwt) + EPS;
            _smw(i,k) = wt;
            wsum += wt;
        }   // end for

        assert( wsum >= EPS);
        _smw.row(i) /= wsum;  // Normalise the row of weights
    }   // end for
}   // end ctor
