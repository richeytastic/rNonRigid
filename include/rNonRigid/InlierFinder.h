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

#ifndef RNONRIGID_INLIER_FINDER_H
#define RNONRIGID_INLIER_FINDER_H

#include "Types.h"

namespace rNonRigid {

class rNonRigid_EXPORT InlierFinder
{
public:
    InlierFinder( float kappa=4.0f, bool useOrientation=true, size_t numIterations=10);

    // Returns a vector giving the probability of each element in flt being an inlier
    // given the updated feature correspondences (crs) and flags.
    VecXf operator()( const MatX6f &flt,      // N rows X 6 columns
                      const MatX6f &crs,      // N rows X 6 columns
                      const VecXf &flags) const; // N rows
private:
    const float _kappa;
    const bool _useOrientation;
    const size_t _numIterations;
    const float _minSig;
    const float _maxSig;
};  // end class

}   // end namespace

#endif
