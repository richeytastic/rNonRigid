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

#include <RigidRegistration.h>
#include <RigidTransformer.h>
using rNonRigid::RigidRegistration;
using rNonRigid::Mat4f;
using rNonRigid::Mesh;

RigidRegistration::RigidRegistration( size_t maxUpdateIts,
                                      size_t k, float flagThresh, bool eqPushPull,
                                      float kappa, bool useOrient, size_t numInlierIts,
                                      bool useScaling)
    :
      _maxUpdateIts( maxUpdateIts),
      _corresponder( k, flagThresh, eqPushPull),
      _inlierFinder( kappa, useOrient, numInlierIts),
      _useScaling( useScaling)
{
}   // end ctor


Mat4f RigidRegistration::operator()( Mesh &flt, const Mesh &tgt, Mat4f nT) const
{
    const K3Tree kdT( tgt.positions());
    const RigidTransformer rgdTrans( _useScaling);

    Mat4f T = Mat4f::Identity();
    VecXf flags;  // Correspondence flags
    for ( size_t i = 0; i < _maxUpdateIts; ++i)
    {
        T = nT * T;
        flt.transform( nT);
        const K3Tree kdF( flt.positions());
        const SparseMat A = _corresponder( kdF, kdT, flags); // F.rows() X T.rows()
        const MatXf crs = A * tgt.features;    // F rows
        const VecXf wts = _inlierFinder( flt.features, crs, flags); // Correspondence weights
        nT = rgdTrans( flt.positions(), crs.leftCols<3>(), wts);  // Calc next transform
        if ( nT.isIdentity( 1e-4f)) // Done if close to not needing another transform
            break;
    }   // end for

    return T;
}   // end operator()
