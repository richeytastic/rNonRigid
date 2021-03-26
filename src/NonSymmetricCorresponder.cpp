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

#include <NonSymmetricCorresponder.h>
using rNonRigid::NonSymmetricCorresponder;
using rNonRigid::SparseMat;
using rNonRigid::K3Tree;


NonSymmetricCorresponder::NonSymmetricCorresponder( size_t k) : _k(k) {}


SparseMat NonSymmetricCorresponder::operator()( const MatX3f &F, const K3Tree& T) const
{
    // knnF2T will iterate over floating and search for correspondences on target
    KNNCorresponder knnF2T( F, _k);
    return normaliseRows( knnF2T.find( T));
}   // end operator()

