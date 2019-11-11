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

#include <KDTree.h>
#include <cassert>
#include <nanoflann.hpp>
using rNonRigid::KDTree;
using rNonRigid::SurfacePoints;
using rNonRigid::FeatMat;
using rNonRigid::FeatVec;
using rNonRigid::NFEATURES;


namespace {
using MyKDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, SurfacePoints<float> >, SurfacePoints<float>, NFEATURES>;
}   // end namespace


class KDTree::Impl
{
public:
    explicit Impl( const FeatMat& m) : _pcloud(m)
    {
        assert( m.cols() == NFEATURES);
        _kdtree = new MyKDTree( m.cols(), _pcloud, nanoflann::KDTreeSingleIndexAdaptorParams(15));
        _kdtree->buildIndex();
    }   // end ctor

    ~Impl() { delete _kdtree;}

    const FeatMat& data() const { return _pcloud.model();}

    size_t findn( const FeatVec& p, size_t n, size_t *nearv, float *sqdis) const
    {
        if ( n == 0)
            return 0;
        return _kdtree->knnSearch( &p(0), n, nearv, sqdis);
    }   // end findn

private:
    const SurfacePoints<float> _pcloud;
    MyKDTree *_kdtree;
};  // end class


KDTree::KDTree( const FeatMat& m) : _impl( new Impl(m)) {}
KDTree::~KDTree() { delete _impl;}
const FeatMat& KDTree::data() const { return _impl->data();}
size_t KDTree::findn( const FeatVec& p, size_t n, size_t *nv, float *sqd) const { return _impl->findn( p, n, nv, sqd);}
