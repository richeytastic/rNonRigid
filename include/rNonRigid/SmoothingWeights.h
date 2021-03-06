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

#ifndef RNONRIGID_SMOOTHING_WEIGHTS_H
#define RNONRIGID_SMOOTHING_WEIGHTS_H

#include "K3Tree.h"

namespace rNonRigid {

class rNonRigid_EXPORT SmoothingWeights
{
public:
    SmoothingWeights( const K3Tree&, size_t K, float sigma);

    const MatXi& indices() const { return _indices;}
    const MatXf& weights() const { return _smw;}

private:
    MatXi _indices;
    MatXf _smw;
};  // end class

}   // end namespace

#endif
