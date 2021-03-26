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

#ifndef RNONRIGID_RIGID_REGISTRATION_H
#define RNONRIGID_RIGID_REGISTRATION_H

#include "SymmetricCorresponder.h"
#include "InlierFinder.h"

namespace rNonRigid {

class rNonRigid_EXPORT RigidRegistration
{
public:
    // Parameters:
    // maxUpdateIts : max iterations before halting regardless of if registration reached.
    // k            : KNN when looking for closest points on the opposite surface.
    // flagThresh   : vertex affinity values higher than this make flag values 1 (all others 0).
    // eqPushPull   : vertex affinities are normalised independently to negate distance bias.
    // kappa        : number of stddevs defining inlier range for discovered correspondences.
    // useOrient    : whether or not to use vertex normals when evaluating inlier correspondences.
    // numInlierIts : number of iterations over which inlier probabilities are re-calculated.
    // useScaling   : transform matrix generated per aligning iteration includes scaling.
    RigidRegistration( size_t maxUpdateIts=200,
                       size_t k=3, float flagThresh=0.9f, bool eqPushPull=true,
                       float kappa=4.0f, bool useOrient=true, size_t numInlierIts=10, 
                       bool useScaling=true);

    // Apply the rigid registration to map mask to target.
    // Optionally provide an initial mask transform.
    // Returns the transform that was applied to mask.
    Mat4f operator()( Mesh &mask, const Mesh &target, Mat4f T=Mat4f::Identity()) const;

private:
    const size_t _maxUpdateIts;
    const SymmetricCorresponder _corresponder;
    const InlierFinder _inlierFinder;
    const bool _useScaling;
};  // end class

}   // end namespace

#endif
