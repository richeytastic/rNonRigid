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

#include <FeaturesCalculate.h>
#include <cassert>
using rNonRigid::FeatMat;
using rNonRigid::FlagVec;
using rNonRigid::SparseMat;


FlagVec rNonRigid::calcFlags( const SparseMat& A, const FlagVec& I, float threshold)
{
    assert( threshold >= 0.0f);
    assert( threshold <= 1.0f);
    FlagVec F = A * I;
    const size_t N = F.size();
    for ( size_t i = 0; i < N; ++i)
        F[i] = F[i] > threshold ? 1.0f : 0.0f;
    return F;
}   // end calcFlags



FeatMat rNonRigid::calcFeatures( const SparseMat& A, const FeatMat& T) { return A * T;}
