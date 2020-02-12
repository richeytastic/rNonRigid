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

#include <RigidRegistration.h>
using rNonRigid::RigidRegistration;
using rNonRigid::Mat4f;
using rNonRigid::Mesh;


RigidRegistration::RigidRegistration( size_t k, float flagThresh, bool eqPushPull,
                                      float kappa, bool useOrient, size_t numInlierIts,
                                      bool useScaling, size_t numUpdateIts)
    : _corresponder( k, flagThresh, eqPushPull),
      _inlierFinder( kappa, useOrient, numInlierIts),
      _transformer( useScaling),
      _numUpdateIts( numUpdateIts)
{
}   // end ctor


Mat4f RigidRegistration::operator()( const Mesh &mask, const Mesh &target) const
{
    Mat4f T = Mat4f::Identity();
    const KDTree kdT( target.features);
    FeatMat F = mask.features;

    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        const KDTree kdF( F);
        FlagVec crsFlags;   // Corresponding target features
        const FeatMat crs = _corresponder( kdF, mask.flags, kdT, target.flags, &crsFlags);

        const VecXf wts = _inlierFinder( F, crs, crsFlags); // Correspondence weightings
        T = _transformer( F, crs, wts) * T; // Update the transform

        // Note transform always from original to reduce error propogation
        F = transformFeatures( mask.features, T);
    }   // end for

    return T;
}   // end operator()
