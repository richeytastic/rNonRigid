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

#include <KNNCorresponder.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
using rNonRigid::KNNCorresponder;
using rNonRigid::KDTree;
using rNonRigid::SparseMat;
using rNonRigid::FeatMat;


void rNonRigid::normaliseRows( SparseMat& m)
{
    const size_t N = m.rows();
    std::vector<float> rowSums( N, 0.0f);

    for ( int i = 0; i < m.outerSize(); ++i)
        for ( SparseMat::InnerIterator it(m,i); it; ++it)
            rowSums[it.row()] += it.value();

    for ( int i = 0; i < m.outerSize(); ++i)
        for ( SparseMat::InnerIterator it(m,i); it; ++it)
            it.valueRef() = it.value() / rowSums[it.row()];
}   // end normaliseRows


KNNCorresponder::KNNCorresponder( const FeatMat& m, size_t k, bool doNorm) : _qry(m), _k(k), _doRowNormalise(doNorm)
{
    assert( k < size_t(m.rows()));
    assert( k >= 1);
}   // end ctor


SparseMat KNNCorresponder::find( const KDTree& kdt) const
{
    const size_t K = _k;
    const size_t n = _qry.rows();          // # query vertices
    const size_t m = kdt.data().rows();    // # target vertices
    assert( K < m);

    static const float EPS = 1e-6f; // Required in the case of any distance == 0 to prevent div-by-zero

    std::vector<size_t> kverts(K);  // K closest vertices on the target model.
    std::vector<float> sqdis(K);    // Corresponding squared distances of each closest vertex to the search vertex
    using Triplet = Eigen::Triplet<float>;
    std::vector<Triplet> aelems( n*K, Triplet(0,0,0.0f));   // Elements for the affinity matrix
    size_t c = 0;   // Counter for aelems

    // For each floating vertex, find the K nearest vertices on the target
    for ( size_t i = 0; i < n; ++i)
    {
        kdt.findn( _qry.row(i), K, &kverts[0], &sqdis[0]);  // Find the k nearest points on the target for query point i
        const Eigen::Vector3f n = _qry.row(i).tail<3>();    // Normal for query point i

        for ( size_t k = 0; k < K; ++k)
        {
            const size_t j = kverts.at(k); // j is the vertex row on the target closest to vertex i of the query set
            float aij = powf( sqdis[k] < EPS ? EPS : sqdis[k], -1); // Affinity weight as inverse squared distance
            // Incorporate the orientation from the matched target vertex
            aij *= 0.5f + n.dot( kdt.data().row(j).tail<3>()) / 2.0f; // Normalise dot product in [0,1]
            // Check for numerical stability since normalizing these elements later and set the entry.
            aelems[c++] = Triplet(i, j, std::max( aij, 1e-4f));
        }   // end for
    }   // end for

    SparseMat A( n, m);
    A.setFromTriplets( aelems.begin(), aelems.end());
    if ( _doRowNormalise)
        normaliseRows( A);

    return A;
}   // end find
