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

#ifndef RNONRIGID_SMOOTHING_WEIGHTS_H
#define RNONRIGID_SMOOTHING_WEIGHTS_H

#include "KNNMap.h"

namespace rNonRigid {

class rNonRigid_EXPORT SmoothingWeights
{
public:
    SmoothingWeights( const KNNMap&, const FlagVec&, float sigma);

    const MatXi& indices() const { return _kmap.indices();}
    const MatXf& weights() const { return _smw;}

private:
    const KNNMap &_kmap;
    MatXf _smw;
};  // end class

}   // end namespace

#endif
