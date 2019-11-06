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

#ifndef RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H
#define RNONRIGID_VISCO_ELASTIC_TRANSFORMER_H

/**
 * This implementation based on MeshMonk registration::ViscoElasticTransformer.
 */

#include "SmoothingWeights.h"

namespace RNonRigid {

class rNonRigid_EXPORT ViscoElasticTransformer
{
public:
    // Parameters:
    ViscoElasticTransformer( size_t numViscousStart, size_t numViscousEnd,
                             size_t numElasticStart, size_t numElasticEnd,
                             size_t numUpdatesTotal,
                             float inlierThresholdWt=0.8f,
                             size_t numOutlierDiffIts=15);

    // Find and return the displacement field iteration to be applied next to the set of features F.
    // Parameters:
    // i    : The iteration. The display field is reset when zero.
    // CF   : Iterative displacement 3-vector field FROM current floating positions TO fresh correspondence points.
    // Fwts : The smoothing weights for the neighbours of the floating vertices in their initial state.
    // iws  : N vector of inlier weights denoting how much each position displacement contributes.
    // flgs : N vector of 1s or 0s denoting which of displacement entries to use.
    // Returns the update displacement vector field.
    DispMat update( size_t i, const DispMat &CF, const SmoothingWeights &Fwts, const VecXf &iws, const FlagVec &flgs);

private:
    const float _viscousAnnealingRate;
    const float _elasticAnnealingRate;
    const size_t _numViscousStart;
    const size_t _numElasticStart;
    const float _inlierThresholdWt;
    const size_t _numOutlierDiffIts;
    DispMat _prevDispField;
    DispMat _thisDispField;
};  // end class

}   // end namespace

#endif
