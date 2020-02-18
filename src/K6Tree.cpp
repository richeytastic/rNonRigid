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

#include <K6Tree.h>
#include <cassert>
#include <nanoflann.hpp>
using rNonRigid::K6Tree;
using rNonRigid::S6Points;
using rNonRigid::FeatMat;
using rNonRigid::FeatVec;


namespace {
using MyK6Tree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, S6Points<float> >, S6Points<float>, 6>;
}   // end namespace


class K6Tree::Impl
{
public:
    explicit Impl( const FeatMat& m) : _pcloud(m)
    {
        assert( m.cols() == 6);
        _kdtree = new MyK6Tree( (int)m.cols(), _pcloud, nanoflann::KDTreeSingleIndexAdaptorParams(15));
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
    const S6Points<float> _pcloud;
    MyK6Tree *_kdtree;
};  // end class


K6Tree::K6Tree( const FeatMat& m) : _impl( new Impl(m)) {}

K6Tree::~K6Tree() { delete _impl;}

const FeatMat& K6Tree::data() const { return _impl->data();}

size_t K6Tree::findn( const FeatVec& p, size_t n, size_t *nv, float *sqd) const
{
    return _impl->findn( p, n, nv, sqd);
}   // end findn
