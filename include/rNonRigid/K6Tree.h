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

#ifndef RNONRIGID_K6_TREE_H
#define RNONRIGID_K6_TREE_H

#include "Types.h"

namespace rNonRigid {

class rNonRigid_EXPORT K6Tree
{
public:
    // Create a K6-tree from the given matrix where the number of rows equate
    // to points and the columns are the coordinates of the vertices.
    explicit K6Tree( const MatX6f&);
    ~K6Tree();

    // Returns the matrix passed in to the constructor.
    const MatX6f& data() const;

    // Returns the number of points in the set.
    inline size_t numPoints() const { return (size_t)data().rows();}

    // Find the n points closest to the given feature. Arrays ridxs and sqdis must
    // be arrays of length n. Returns actual number of points found which may be less than n.
    size_t findn( const FeatVec&, size_t n, size_t *ridxs, float *sqdis) const;

private:
    class Impl;
    Impl *_impl;

    K6Tree( const K6Tree&) = delete;
    K6Tree& operator=( const K6Tree&) = delete;
};  // end class


template <typename T>
class S6Points
{
public:
    S6Points( const MatX6f& m) : _m(m) {}
    inline size_t kdtree_get_point_count() const { return (size_t)_m.rows();}
    inline T kdtree_get_pt( const size_t idx, int dim) const { return _m((int)idx, dim);}
    template <class BBOX>
    inline bool kdtree_get_bbox( BBOX&) const { return false;}
    const MatX6f& model() const { return _m;}
private:
    const MatX6f _m;
};  // end class

}   // end namespace

#endif
