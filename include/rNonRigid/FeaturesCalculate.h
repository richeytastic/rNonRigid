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

#ifndef RNONRIGID_FEATURES_CALCULATE_H
#define RNONRIGID_FEATURES_CALCULATE_H

#include "Types.h"

namespace rNonRigid {

// Calculates the weighted sum over the target points for new query points and returns as new features set with Q rows and X columns.
// Affine matrix A has Q rows and T columns (Query points, Target points) and must be row normalised.
// Matrix T is the set of target features having T rows and X columns.
rNonRigid_EXPORT FeatMat calcFeatures( const SparseMat& A, const FeatMat& T);

// Creates and returns a new flags vector with Q entries (with entries either 1.0 or 0.0).
// Affinity matrix as before; flags is a column vector with T entries specifying which target points to consider and should
// be set with either 1s or 0s. In the returned vector, all values will be zero except values > threshold which are 1.0.
rNonRigid_EXPORT FlagVec calcFlags( const SparseMat& A, const FlagVec& flags, float threshold=0.5f);

}   // end namespace

#endif
