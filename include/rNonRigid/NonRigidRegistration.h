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

#ifndef RNONRIGID_NON_RIGID_REGISTRATION_H
#define RNONRIGID_NON_RIGID_REGISTRATION_H

/**
 * This implementation based on MeshMonk registration::NonRigidRegistration.
 */

#include "SymmetricCorresponder.h"
#include "InlierFinder.h"
#include "ViscoElasticTransformer.h"

namespace rNonRigid {

class rNonRigid_EXPORT NonRigidRegistration
{
public:
    // Parameters:
    // k                : KNN when looking for closest points on the opposite surface.
    // flagThresh       : affinity values higher than this make flag values 1 (all others 0).
    // eqPushPull       : symmetric correspondence affinity matrices independently row normalised before merging.
    // kappa            : number of stddevs defining inlier range for discovered correspondences.
    // useOrient        : whether or not to use vertex normals when evaluating inlier correspondences.
    // numInlierIts     : number of iterations over which inlier probabilities are re-calculated.
    // smoothK          : number of vertex neighbours in the local smoothing region for each floating vertex.
    // sigmaSmooth      : smoothing Gaussian width for the floating points.
    // numViscousStart  : starting number of viscous steps when beginning transform.
    // numViscousEnd    : final number of viscous steps when finishing transform.
    // numElasticStart  : starting number of elastic steps when beginning transform.
    // numElasticEnd    : final number of elastic steps when finishing transform.
    // numUpdateIts     : number of iterations of the update loop for recalculating/applying correspondences.
    NonRigidRegistration( size_t k=3, float flagThresh=0.9f, bool eqPushPull=false,
                          float kappa=4.0f, bool useOrient=true, size_t numInlierIts=10,
                          size_t smoothK=80, float sigmaSmooth=3.0f,
                          size_t numViscousStart=100, size_t numViscousEnd=1,
                          size_t numElasticStart=100, size_t numElasticEnd=1,
                          size_t numUpdateIts=200);

    // Find the non-rigid registration between F and T where points are stored row wise
    // with each row having 6 elements as X,Y,Z position and X,Y,Z normal.
    // Parameters:
    // mask   : the floating mesh - updated to mapped position on return.
    // target : the target mesh to which the mask will be registered.
    void operator()( Mesh &mask, const Mesh &target) const;

private:
    const size_t _smoothK;
    const float _sigmaSmooth;
    const size_t _numUpdateIts;
    const SymmetricCorresponder _corresponder;
    const InlierFinder _inlierFinder;
    mutable ViscoElasticTransformer _transformer;
};  // end class

}   // end namespace

#endif
