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

#ifndef RNONRIGID_KNN_CORRESPONDER_H
#define RNONRIGID_KNN_CORRESPONDER_H

/**
 * Does distance weighted K nearest neighbour to find the nearest points on a target mesh
 * corresponding to a given query mesh. Returns an affinity matrix of correspondences
 * with values in proportion to the inverse of the squared distances between the points.
 */
#include "KDTree.h"

namespace RNonRigid {

class rNonRigid_EXPORT KNNCorresponder
{
public:
    // Set the query points as the first three columns of the given matrix.
    // k is the number of nearest neighbours on the target to search for.
    KNNCorresponder( const FeatMat& query, size_t k=3, bool doRowNormalise=true);

    // Return the Q x T affinity matrix where Q is the number of points in the query set
    // and T the number of points in the target set. Each entry is the inverse of the squared
    // distance between the query and the target points. Only vertices from the target set
    // found within K nearest neighbours of the query points have non-zero entries.
    // Correspondence points C can be calculated from provided target points T and the
    // returned matrix A as C = AT. That is, each point coregistered to the query set
    // is the weighted sum of points in the target set.
    SparseMat find( const KDTree& target) const;

private:
    const FeatMat& _qry;
    const size_t _k;
    const bool _doRowNormalise;
};  // end class


// Normalise the rows of the given sparse matrix.
rNonRigid_EXPORT void normaliseRows( SparseMat&);

}   // end namespace

#endif
