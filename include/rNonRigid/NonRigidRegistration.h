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

#ifndef RNONRIGID_NON_RIGID_REGISTRATION_H
#define RNONRIGID_NON_RIGID_REGISTRATION_H

#include "SymmetricCorresponder.h"
#include "InlierFinder.h"

namespace rNonRigid {

class rNonRigid_EXPORT NonRigidRegistration
{
public:
    // Parameters:
    // numUpdateIts     : number of iterations required.
    // k                : KNN when looking for closest points on the opposite surface.
    // flagThresh       : affinity values higher than this make flag values 1 (all others 0).
    // eqPushPull       : vertex affinities are normalised independently to negate distance bias.
    // kappa            : number of stddevs defining inlier range for discovered correspondences.
    // useOrient        : whether or not to use vertex normals when evaluating inlier correspondences.
    // numInlierIts     : # iterations over which inlier probabilities are re-calculated.
    // smoothK          : # vertex neighbours in local smoothing region for each floating vertex.
    // smoothS          : smoothing Gaussian width for the floating points.
    // numViscousStart  : starting number of viscous steps when beginning transform.
    // numViscousEnd    : final number of viscous steps when finishing transform.
    // numElasticStart  : starting number of elastic steps when beginning transform.
    // numElasticEnd    : final number of elastic steps when finishing transform.
    NonRigidRegistration( size_t numUpdateIts=200,
                          size_t k=3, float flagThresh=0.9f, bool eqPushPull=true,
                          float kappa=4.0f, bool useOrient=true, size_t numInlierIts=10,
                          size_t smoothK=80, float smoothS=3.0f,
                          size_t numViscousStart=100, size_t numViscousEnd=1,
                          size_t numElasticStart=100, size_t numElasticEnd=1);

    // Find the non-rigid registration between F and T where points are stored row
    // wise with each row having 6 elements as X,Y,Z position and X,Y,Z normal.
    // On return, F has its features registered to surface T.
    // F : the floating template to map to the target.
    // T : the target to which the floating template is mapped.
    void operator()( Mesh &F, const Mesh &T) const;

private:
    const size_t _numUpdateIts;
    const size_t _smoothK;
    const float _smoothS;
    const SymmetricCorresponder _corresponder;
    const InlierFinder _inlierFinder;
    const size_t _nvStart, _nvEnd;
    const size_t _neStart, _neEnd;
};  // end class

}   // end namespace

#endif
