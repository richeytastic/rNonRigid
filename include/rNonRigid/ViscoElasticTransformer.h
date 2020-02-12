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

#ifndef RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H
#define RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H

/**
 * This implementation based on MeshMonk registration::ViscoElasticTransformer.
 */

#include "SmoothingWeights.h"

namespace rNonRigid {

class rNonRigid_EXPORT ViscoElasticTransformer
{
public:
    // Parameters:
    ViscoElasticTransformer( size_t numViscousStart, size_t numViscousEnd,
                             size_t numElasticStart, size_t numElasticEnd,
                             size_t numUpdatesTotal,
                             float inlierThresholdWt=0.8f,
                             size_t numOutlierDiffIts=15);

    // Modify the given displacement field.
    // Parameters:
    // i    : The update iteration.
    // D    : Displacement field to update.
    // swts : The smoothing weights for the neighbours of the floating vertices in their initial state.
    // iwts : N vector of inlier weights denoting how much each position displacement contributes.
    // Returns the updated displacement vector field.
    void update( size_t i, MatX3f &D, const SmoothingWeights &swts, const VecXf &iwts);

private:
    const float _viscousAnnealingRate;
    const float _elasticAnnealingRate;
    const size_t _numViscousStart;
    const size_t _numElasticStart;
    const float _inlierThresholdWt;
    const size_t _numOutlierDiffIts;
    MatX3f _prevField;
    MatX3f _thisField;
};  // end class


rNonRigid_EXPORT MatXf createWeights( const SmoothingWeights&, const VecXf&);

rNonRigid_EXPORT void regularise( MatX3f&, const MatXf&, const MatXi&, size_t);

rNonRigid_EXPORT void diffuseOutliers( MatX3f&, const SmoothingWeights&, const VecXf&, float, size_t);

}   // end namespace

#endif
