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

#include <SymmetricCorresponder.h>
#include <cassert>
using rNonRigid::SymmetricCorresponder;
using rNonRigid::SparseMat;
using rNonRigid::K3Tree;
using rNonRigid::VecXf;


SymmetricCorresponder::SymmetricCorresponder( size_t k, float h, bool eqpp)
    : _k(k), _thresh(h), _eqpp(eqpp)
{
    assert( h >= 0.0f);
    assert( h <= 1.0f);
}   // end ctor


namespace {
VecXf calcFlags( const SparseMat &A, const VecXf &f, float threshold)
{
    assert( A.cols() == f.rows());
    VecXf F = A * f;
    for ( long i = 0; i < F.size(); ++i)
        F[i] = F[i] > threshold ? 1.0f : 0.0f;
    assert( F.size() == A.rows());
    return F;
}   // end calcFlags
}   // end namespace


SparseMat SymmetricCorresponder::operator()( const K3Tree& F, const K3Tree& T, VecXf &fC) const
{
    // knnF2T will iterate over floating and search for correspondences on target
    // knnT2F will iterate over target and search for correspondences on floating
    const KNNCorresponder knnF2T( F.data(), _k);  // Push floating to target
    const KNNCorresponder knnT2F( T.data(), _k);  // Pull floating to target

    // For F vertices in the floating set, and T vertices in the target set
    const SparseMat A_ft = knnF2T.find( T);   // Affinity matrix F x T (not row normalised)
    const SparseMat A_tf = knnT2F.find( F);   // Affinity matrix T x F (not row normalised)
    const SparseMat A_ft_n = normaliseRows(A_ft);
    const SparseMat A_tf_n = normaliseRows(A_tf);

    SparseMat A;
    if ( _eqpp)
    {
        A = normaliseRows( A_ft_n + SparseMat( A_tf_n.transpose()));
        fC = VecXf::Ones(A.cols());
    }   // end if
    else
    {
        // Check distance thresholds on each lookup independently (requires normalised)
        fC = calcFlags( A_ft_n, VecXf::Ones(A_ft_n.cols()), _thresh);
        fC = calcFlags( A_tf_n, fC, _thresh);
        A = normaliseRows( A_ft + SparseMat( A_tf.transpose()));
        assert( fC.size() == F.data().rows());
    }   // end else

    fC = calcFlags( A, fC, _thresh);

    return A;
}   // end operator()

