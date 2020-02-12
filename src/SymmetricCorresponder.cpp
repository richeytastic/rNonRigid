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

#include <SymmetricCorresponder.h>
#include <cassert>
using rNonRigid::SymmetricCorresponder;
using rNonRigid::KDTree;
using rNonRigid::SparseMat;
using rNonRigid::FeatMat;
using rNonRigid::FlagVec;

namespace {
FlagVec calcFlags( const SparseMat& A, const FlagVec& I, float threshold)
{
    FlagVec F = A * I;
    const size_t N = F.size();
    for ( size_t i = 0; i < N; ++i)
        F[i] = F[i] > threshold ? 1.0f : 0.0f;
    return F;
}   // end calcFlags
}   // end namespace


SymmetricCorresponder::SymmetricCorresponder( size_t k, float h, bool eqpp)
    : _k(k), _thresh(h), _eqpp(eqpp)
{
    assert( h >= 0.0f);
    assert( h <= 1.0f);
}   // end ctor


FeatMat SymmetricCorresponder::operator()( const KDTree& F, const FlagVec& fF,
                                           const KDTree& T, const FlagVec& fT,
                                           FlagVec *fC) const
{
    //assert( fF.size() == (int)F.numPoints());
    //assert( fT.size() == (int)T.numPoints());

    // knnF2T will iterate over floating and search for correspondences on target
    // knnT2F will iterate over target and search for correspondences on floating
    KNNCorresponder knnF2T( F.data(), _k, false);  // Push floating to target
    KNNCorresponder knnT2F( T.data(), _k, false);  // Pull floating to target

    // For F vertices in the floating set, and T vertices in the target set
    SparseMat A = knnF2T.find( T);   // Affinity matrix F x T (not row normalised)
    SparseMat B = knnT2F.find( F);   // Affinity matrix T x F (not row normalised)

    // A and B are not row normalised, but normalised versions
    // are needed for setting the per affinity matrix flags.
    SparseMat Anorm, Bnorm;
    if ( fC || _eqpp)
    {
        Anorm = A;
        Bnorm = B;
        normaliseRows( Anorm);
        normaliseRows( Bnorm);
    }   // end if

    FlagVec flags;
    if ( fC)
    {
        // Get flags per individual affinity matrix using input flags of both floating and target
        // sets to mask which of the common points to use for generating new correspondences.
        flags = calcFlags( Anorm, fT, _thresh); // flags.size() == A.rows()
        flags = calcFlags( Bnorm, flags.cwiseProduct( fF), _thresh); // flags.size() == A.cols()
    }   // end if

    // If equalising matrices then magnitudes are irrelevant so
    // set A and B to be the normalised versions prior to merging.
    if ( _eqpp)
    {
        A = Anorm;
        B = Bnorm;
    }   // end if

    A += SparseMat( B.transpose()); // Merge
    normaliseRows(A);

    if ( fC)    // Update flags?
        *fC = calcFlags( A, flags, _thresh);

    return A * T.data();    // Return corresponding features
}   // end operator()

