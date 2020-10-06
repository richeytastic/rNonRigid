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

#ifndef RNONRIGID_FAST_DEFORM_REGISTRATION_H
#define RNONRIGID_FAST_DEFORM_REGISTRATION_H

#include "Types.h"

namespace rNonRigid {

class rNonRigid_EXPORT FastDeformRegistration
{
public:
    // Parameters:
    // n                : Number of iterations
    explicit FastDeformRegistration( size_t n);

    // Deform mask to target.
    void operator()( Mesh &mask, const Mesh &target) const;

private:
    const size_t _numUpdateIts;
};  // end class

}   // end namespace

#endif
