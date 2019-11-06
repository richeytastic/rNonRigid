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

#include <RigidRegistration.h>
#include <iostream>
#include <ctime>
using RNonRigid::RigidRegistration;
using RNonRigid::Mat4f;
using RNonRigid::FeatMat;
using RNonRigid::FlagVec;

RigidRegistration::RigidRegistration( size_t k, float flagThresh, bool eqPushPull,
                                      float kappa, bool useOrient, size_t numInlierIts,
                                      bool useScaling, size_t numUpdateIts)
    : _corresponder( k, flagThresh, eqPushPull),
      _inlierFinder( kappa, useOrient, numInlierIts),
      _transformer( useScaling),
      _numUpdateIts( numUpdateIts)
{
}   // end ctor


Mat4f RigidRegistration::operator()( FeatMat &uF, const FlagVec &fFlags, const FeatMat &tgt, const FlagVec &tFlags) const
{
    Mat4f T = Mat4f::Identity();
    const KDTree tkdt( tgt);
    const FeatMat F = uF; // Const version of F for reliable (non-error propogating) iterative application of transform

    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        FlagVec cFlags;
        const FeatMat crs = _corresponder( uF, fFlags, tkdt, tFlags, cFlags); // Returns corresponding features on target
        const VecXf wts = _inlierFinder( uF, crs, cFlags);          // Returns correspondence weightings
        T = _transformer( uF, crs, wts) * T;                        // Update the transform
        uF = transformFeatures( F, T);  // Note transform always from original F to avoid propogation of errors
    }   // end for

    return T;
}   // end operator()
