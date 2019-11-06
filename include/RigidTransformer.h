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

#ifndef RNONRIGID_RIGID_TRANSFORMER_H
#define RNONRIGID_RIGID_TRANSFORMER_H

/**
 * This implementation based on MeshMonk registration::RigidTransformer.
 * Note that unlike the MeshMonk implementation, RigidTransformer::find
 * does not apply the transform to the input points (use helper function
 * transformFeatures for that).
 */

#include "Types.h"

namespace RNonRigid {

class rNonRigid_EXPORT RigidTransformer
{
public:
    // Parameters:
    // useScaling   : the transformation matrix returned from find will incorporate a scaling estimate.
    explicit RigidTransformer( bool useScaling=true);

    // Find and return the rigid transform between two sets of features with weights.
    // Parameters:
    // F  : N rows with vertex X,Y,Z positions in the first three columns - other columns are ignored.
    // C  : N rows with vertex X,Y,Z positions in the first three columns - other columns are ignored.
    // w  : vector of floats denoting how much each feature mapping contributes to the transform.
    // Returns the transform matrix M to get to C from F.
    Mat4f operator()( const FeatMat &F, const FeatMat &C, const VecXf &w) const;

private:
    const bool _useScaling;
};  // end class


// Apply T to the features matrix F which should be 6 elements long (columns) in each row with the
// first three elements giving a point's position, and the last three elements giving its normal.
rNonRigid_EXPORT FeatMat transformFeatures( const FeatMat& F, const Mat4f& T);

}   // end namespace

#endif
