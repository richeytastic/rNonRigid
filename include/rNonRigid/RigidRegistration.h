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

#ifndef RNONRIGID_RIGID_REGISTRATION_H
#define RNONRIGID_RIGID_REGISTRATION_H

/**
 * This implementation based on MeshMonk registration::RigidRegistration.
 * Non-symmetric correspondence is not implemented.
 */

#include "SymmetricCorresponder.h"
#include "InlierFinder.h"
#include "RigidTransformer.h"

namespace rNonRigid {

class rNonRigid_EXPORT RigidRegistration
{
public:
    // Parameters:
    // k            : KNN when looking for closest points on the opposite surface.
    // flagThresh   : affinity values higher than this make flag values 1 (all others 0).
    // eqPushPull   : symmetric correspondence affinity matrices independently row normalised before merging.
    // kappa        : number of stddevs defining inlier range for discovered correspondences.
    // useOrient    : whether or not to use vertex normals when evaluating inlier correspondences.
    // numInlierIts : number of iterations over which inlier probabilities are re-calculated.
    // useScaling   : transform matrix generated per aligning iteration includes scaling.
    // numUpdateIts : number of iterations of the update loop for recalculating/applying correspondences.
    RigidRegistration( size_t k=3, float flagThresh=0.9f, bool eqPushPull=false,
                       float kappa=4.0f, bool useOrient=true, size_t numInlierIts=10, 
                       bool useScaling=true, size_t numUpdateIts=20);

    // Find the rigid registration between mask and target.
    // Returns the transformation matrix that should be applied to the mask.
    // Parameters:
    // mask     : the floating mesh for which a transform to the target will be found.
    // target   : the target mesh to which the mask will be registered.
    Mat4f operator()( const Mesh &mask, const Mesh &target) const;

private:
    const SymmetricCorresponder _corresponder;
    const InlierFinder _inlierFinder;
    const RigidTransformer _transformer;
    const size_t _numUpdateIts;
};  // end class

}   // end namespace

#endif
