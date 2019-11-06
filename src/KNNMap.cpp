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

#include <KNNMap.h>
#include <KDTree.h>
using RNonRigid::KNNMap;
using RNonRigid::FeatMat;


KNNMap::KNNMap( const FeatMat &qry, const FeatMat &tgt, size_t K)
    : _idxs( qry.rows(), K), _sqds( qry.rows(), K)
{
    const KDTree kdt( tgt);
    const size_t N = qry.rows();
    std::vector<size_t> idxs(K);
    std::vector<float> sqds(K);
    for ( size_t i = 0; i < N; ++i)
    {
        kdt.findn( qry.row(i), K, &idxs[0], &sqds[0]);
        for ( size_t k = 0; k < K; ++k)
        {
            _idxs(i,k) = int(idxs[k]);
            _sqds(i,k) = sqds[k];
        }   // end for
    }   // end for
}   // end ctor


