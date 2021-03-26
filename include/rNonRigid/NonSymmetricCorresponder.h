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

#ifndef RNONRIGID_NON_SYMMETRIC_CORRESPONDER_H
#define RNONRIGID_NON_SYMMETRIC_CORRESPONDER_H

#include "KNNCorresponder.h"

namespace rNonRigid {

class rNonRigid_EXPORT NonSymmetricCorresponder
{
public:
    // k : each point in F looks for k nearest neighbours on T
    explicit NonSymmetricCorresponder( size_t k=3);

    // Find and return affinity matrix A between F and T. Used to calculate a set
    // of features as A * T.data() corresponding to the entries of F.
    // F      : F rows by 3 columns query points (floating mask)
    // T      : kd-tree for the target points (has T points)
    SparseMat operator()( const MatX3f& F, const K3Tree& T) const;

private:
    size_t _k;
};  // end class

}   // end namespace

#endif


