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

#include <NonRigidRegistration.h>
#include <ViscoElasticTransformer.h>
#include <SmoothingWeights.h>
#include <memory>
using rNonRigid::NonRigidRegistration;
using rNonRigid::Mesh;


NonRigidRegistration::NonRigidRegistration( size_t numUpdateIts,
                                            size_t k, float flagThresh, bool eqPushPull,
                                            float kappa, bool useOrient, size_t numInlierIts,
                                            size_t smoothK, float smoothS,
                                            size_t nvStart, size_t nvEnd,
                                            size_t neStart, size_t neEnd)
    :
      _numUpdateIts( numUpdateIts),
      _smoothK( smoothK), _smoothS( smoothS),
      _corresponder( k, flagThresh, eqPushPull),
      _inlierFinder( kappa, useOrient, numInlierIts),
      _nvStart(nvStart), _nvEnd(nvEnd),
      _neStart(neStart), _neEnd(neEnd)
{
}   // end ctor


void NonRigidRegistration::operator()( Mesh &flt, const Mesh &tgt) const
{
    // The KD-tree for the floating surface is rebuilt every iteration...
    std::shared_ptr<K3Tree> kdF = std::shared_ptr<K3Tree>( new K3Tree( flt.features.leftCols(3)));
    const K3Tree kdT( tgt.features.leftCols(3));  // ...while the target is unchanging.

    // Only need to define the smoothing weights once for the floating surface since each vertex
    // weight is calculated only from the distribution of locally neighbouring points and this
    // isn't based on distance (which is updated with each iteration).
    const SmoothingWeights smw( *kdF, _smoothK, _smoothS);

    ViscoElasticTransformer vetrans( smw, _nvStart, _nvEnd, _neStart, _neEnd, _numUpdateIts);

    VecXf flags;  // Correspondence flags updated every iteration by the symmetric corresponder
    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        const SparseMat A = _corresponder( *kdF, kdT, flags);   // F.rows() X T.rows()
        assert( flags.size() == flt.features.rows());
        const MatX6f crs = A * tgt.features;   // F rows
        const VecXf wts = _inlierFinder( flt.features, crs, flags); // Correspondence weights

        // Displacement field from current mask points to corresponding points on tgt
        MatX3f df = crs.leftCols(3) - flt.features.leftCols(3);
        vetrans.update( df, wts); // Regularise, add, then relax back total deformation field.
        flt.update( df);    // Update

        if ( i < _numUpdateIts - 1)
            kdF = std::shared_ptr<K3Tree>( new K3Tree( flt.features.leftCols(3)));
    }   // end for
}   // end operator()
