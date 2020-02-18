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

#include <K3Tree.h>
#include <cassert>
#include <nanoflann.hpp>
using rNonRigid::K3Tree;
using rNonRigid::S3Points;
using rNonRigid::MatX3f;
using rNonRigid::Vec3f;


namespace {
using MyK3Tree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, S3Points<float> >, S3Points<float>, 3>;
}   // end namespace


class K3Tree::Impl
{
public:
    explicit Impl( const MatX3f &m) : _pcloud(m)
    {
        assert( m.cols() == 3);
        _kdtree = new MyK3Tree( (int)m.cols(), _pcloud, nanoflann::KDTreeSingleIndexAdaptorParams(15));
        _kdtree->buildIndex();
    }   // end ctor

    ~Impl() { delete _kdtree;}

    const MatX3f& data() const { return _pcloud.model();}

    size_t findn( const Vec3f& p, size_t n, size_t *nearv, float *sqdis) const
    {
        if ( n == 0)
            return 0;
        return _kdtree->knnSearch( &p(0), n, nearv, sqdis);
    }   // end findn

private:
    const S3Points<float> _pcloud;
    MyK3Tree *_kdtree;
};  // end class


K3Tree::K3Tree( const MatX3f& m) : _impl( new Impl(m)) {}

K3Tree::~K3Tree() { delete _impl;}

const MatX3f& K3Tree::data() const { return _impl->data();}

size_t K3Tree::findn( const Vec3f& p, size_t n, size_t *nv, float *sqd) const
{
    return _impl->findn( p, n, nv, sqd);
}   // end findn
