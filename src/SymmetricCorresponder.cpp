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

#include <SymmetricCorresponder.h>
#include <FeaturesCalculate.h>
#include <cassert>
using RNonRigid::SymmetricCorresponder;
using RNonRigid::KDTree;
using RNonRigid::SparseMat;
using RNonRigid::FeatMat;
using RNonRigid::FlagVec;

SymmetricCorresponder::SymmetricCorresponder( size_t k, float h, bool eqpp) : _k(k), _thresh(h), _eqpp(eqpp) {}


FeatMat SymmetricCorresponder::operator()( const FeatMat& qdata, const FlagVec& qflags, const KDTree& tkdt, const FlagVec& tflags, FlagVec& cflags) const
{
    assert( qflags.size() == qdata.rows());
    assert( tflags.size() == (int)tkdt.numPoints());

    // knnQ2T iterates over the query points and searches for correspondences on the target points
    // knnT2Q iterates over the target points and searches for correspondences on the query points
    KNNCorresponder knnQ2T( qdata, _k, false);         // Push query to target
    KNNCorresponder knnT2Q( tkdt.data(), _k, false);   // Pull query to target

    // For Q vertices in the query set, and T vertices in the target set
    SparseMat A = knnQ2T.find( tkdt);   // Affinity matrix Q x T (not row normalised)
    const KDTree qkdt( qdata);
    SparseMat B = knnT2Q.find( qkdt);   // Affinity matrix T x Q (not row normalised)

    assert( A.rows() == B.cols());
    assert( A.cols() == B.rows());

    // A and B are not row normalised, but normalised versions are needed for setting the per affinity matrix flags.
    SparseMat Anorm = A;
    SparseMat Bnorm = B;
    normaliseRows( Anorm);
    normaliseRows( Bnorm);

    // Get the flags per individual affinity matrix using the input flags to both the query and target
    // sets to mask which of the common points we will use for generating new correspondences.
    FlagVec flags = calcFlags( Anorm, tflags);
    assert( flags.size() == A.rows());
    flags = calcFlags( Bnorm, flags.cwiseProduct( qflags));
    assert( flags.size() == A.cols());

    // If equalising the matrices then magnitudes are irrelevant so we set A and B to be the normalised versions prior to merging.
    if ( _eqpp)
    {
        A = Anorm;
        B = Bnorm;
    }   // end if

    A += SparseMat( B.transpose()); // Merge
    normaliseRows(A);

    // Update flags and return corresponding features
    cflags = calcFlags( A, flags, _thresh);
    return calcFeatures( A, tkdt.data());
}   // end operator()

