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

#include <SmoothingWeights.h>
#include <KNNMap.h>
#include <cfloat>
#include <cmath>
using rNonRigid::SmoothingWeights;
using rNonRigid::K3Tree;


SmoothingWeights::SmoothingWeights( const K3Tree &kdt, size_t K, float sigma)
{
    const size_t N = kdt.data().rows();
    const KNNMap kmap( kdt.data(), kdt, K);
    _indices = kmap.indices();
    _smw = MatXf( N,K);

    static const float EPS = FLT_MIN;
    static const float ONE_MINUS_EPS = 1.0f - EPS;
    const float EXP_FACTOR = -0.5f / std::pow( sigma, 2.0f);

    for ( size_t i = 0; i < N; ++i)
    {
        float wsum = 0.0f;
        for ( size_t k = 0; k < K; ++k)
        {
            const float dsq = kmap.sqDiffs()(i,k);      // Squared distance
            const float gwt = std::exp( dsq * EXP_FACTOR);  // Gaussian weight
            // Rescale weight in [EPS,1] instead of [0,1] so that nodes with an inlier weight
            // of 0 don't end up with a deformation vector having zero magnitude.
            const float wt = ONE_MINUS_EPS * gwt + EPS;
            _smw(i,k) = wt;
            wsum += wt;
        }   // end for
        _smw.row(i) /= wsum;  // Normalise the row of weights
    }   // end for
}   // end ctor
