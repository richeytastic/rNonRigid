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

#ifndef RNONRIGID_SYMMETRIC_CORRESPONDER_H
#define RNONRIGID_SYMMETRIC_CORRESPONDER_H

#include "KNNCorresponder.h"

namespace rNonRigid {

class rNonRigid_EXPORT SymmetricCorresponder
{
public:
    // k          : each point in Q looks for k nearest neighbours on T and vice versa
    // flagThresh : affinity values higher than this cause flag values to be 1 (all others 0).
    // eqPushPull : push and pull affinity matrices are first independently row normalised before merging
    //              so that calculated features and flags do not bias either the query or target points.
    SymmetricCorresponder( size_t k=3, float flagThresh=0.9f, bool eqPushPull=false);

    // Calculate the affinity matrix A between Q and T and use this to calculate and return a set
    // of features on T that correspond to the entries (rows) of Q. Also sets flags for these
    // corresponding features as out parameter cflags.
    // Q          : matrix with Q rows 6 columns.
    // Qflags     : Q vector with elements in {0,1}; only elements with 1 are used to update correspondences.
    // T          : kd-tree for the target points (has T points)
    // Tflags     : T vector with elements in {0,1}; only elements with 1 are used to update correspondences.
    FeatMat operator()( const FeatMat& Q, const FlagVec& Qflags, const KDTree& T, const FlagVec& Tflags, FlagVec& cflags) const;

private:
    size_t _k;
    float _thresh;
    bool _eqpp;
};  // end class

}   // end namespace

#endif
