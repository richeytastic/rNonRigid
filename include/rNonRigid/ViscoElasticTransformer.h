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

#ifndef RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H
#define RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H

#include "SmoothingWeights.h"

namespace rNonRigid {

class rNonRigid_EXPORT ViscoElasticTransformer
{
public:
    // swts : The smoothing weights for the neighbours of the floating vertices in their initial state.
    ViscoElasticTransformer( const SmoothingWeights &swts,
                             size_t numViscousStart, size_t numViscousEnd,
                             size_t numElasticStart, size_t numElasticEnd,
                             size_t numUpdatesTotal,
                             float inlierThresholdWt=0.8f,
                             size_t numOutlierDiffIts=15);

    // Update the displacement field to add for the iteration.
    // iwts : N vector of inlier weights denoting how much each displacement contributes.
    void update( MatX3f&, const VecXf &iwts);

private:
    const SmoothingWeights &_swts;
    const float _viscousAnnealingRate;
    const float _elasticAnnealingRate;
    const size_t _numViscousStart;
    const size_t _numElasticStart;
    const float _inlierThresholdWt;
    const size_t _numOutlierDiffIts;
    MatX3f _field;
    float _i;
};  // end class

}   // end namespace

#endif
