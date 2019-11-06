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

#ifndef RNONRIGID_KNN_MAP_H
#define RNONRIGID_KNN_MAP_H

#ifdef _WIN32
#pragma warning( disable : 4251)
#endif

/**
 * Find the nearest k points (rows) on the target matrix for every point (row)
 * in the query matrix. The target and query matrices can be the same.
 */
#include "Types.h"

namespace RNonRigid {

class rNonRigid_EXPORT KNNMap
{
public:
    KNNMap( const FeatMat &query, const FeatMat &target, size_t k);

    const MatXi& indices() const { return _idxs;}
    const MatXf& sqDiffs() const { return _sqds;}

private:
    MatXi _idxs;
    MatXf _sqds;
};  // end class

}   // end namespace

#endif
