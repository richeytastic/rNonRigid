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

#ifndef RNONRIGID_RIGID_TRANSFORMER_H
#define RNONRIGID_RIGID_TRANSFORMER_H

#include "Types.h"

namespace rNonRigid {

class rNonRigid_EXPORT RigidTransformer
{
public:
    // useScaling   : the transformation matrix returned from find will incorporate a scaling estimate.
    explicit RigidTransformer( bool useScaling=true);

    // Find and return the rigid transform between two sets of features with weights.
    // F  : N rows with vertex X,Y,Z positions in the columns.
    // C  : N rows with vertex X,Y,Z positions in the columns.
    // w  : N floats in [0,1] denoting how much each position contributes to the transform.
    // Returns transform matrix M to get to C from F.
    Mat4f operator()( const MatX3f &F, const MatX3f &C, const VecXf &w) const;

private:
    const bool _useScaling;
};  // end class


}   // end namespace

#endif
