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
using rNonRigid::MatX6f;
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


namespace {
using namespace rNonRigid;
void _transformFeatures( MatX6f& F, const Mat4f& T)
{
    MatX4f fps( F.rows(), 4);
    fps.leftCols(3) = F.leftCols(3);
    fps.rightCols(1) = VecXf::Ones( F.rows());
    F.leftCols(3) = (fps * T.transpose()).leftCols(3);
    const Mat3f R = T.block<3,3>(0,0) / T(3,3);  // Rotation submatrix
    F.rightCols(3) = F.rightCols(3) * R.transpose();
}   // end _transformFeatures
}   // end namespace


rNonRigid::Mat4f RigidRegistration::operator()( Mesh &flt, const Mesh &tgt, Mat4f T) const
{
    Mat4f nT = T;
    T = Mat4f::Identity();

    const K3Tree kdT( tgt.features.leftCols(3));
    const RigidTransformer rgdTrans( _useScaling);

    VecXf flags;  // Correspondence flags
    for ( size_t i = 0; i < _maxUpdateIts; ++i)
    {
        T = nT * T;
        _transformFeatures( flt.features, nT);
        const K3Tree kdF( flt.features.leftCols(3));

        const SparseMat A = _corresponder( kdF, kdT, flags); // F.rows() X T.rows()
        const MatX6f crs = A * tgt.features;    // F rows
        const VecXf wts = _inlierFinder( flt.features, crs, flags); // Correspondence weights
        nT = rgdTrans( flt.features.leftCols(3), crs.leftCols(3), wts);  // Calc next transform
        if ( nT.isIdentity( 1e-4f)) // Done if close to not needing another transform
            break;
    }   // end for

    return T;
}   // end operator()
