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

#include <NonRigidRegistration.h>
#include <SmoothingWeights.h>
#include <iostream>
#include <memory>
using rNonRigid::NonRigidRegistration;
using rNonRigid::Mesh;


NonRigidRegistration::NonRigidRegistration( size_t k, float flagThresh, bool eqPushPull,
                                            float kappa, bool useOrient, size_t numInlierIts,
                                            size_t smoothK, float sigmaS,
                                            size_t nvStart, size_t nvEnd,
                                            size_t neStart, size_t neEnd,
                                            size_t numUpdateIts)
    : _smoothK( smoothK),
      _sigmaSmooth( sigmaS),
      _numUpdateIts( numUpdateIts),
      _corresponder( k, flagThresh, eqPushPull),
      _inlierFinder( kappa, useOrient, numInlierIts),
      _transformer( nvStart, nvEnd, neStart, neEnd, numUpdateIts)
{
}   // end ctor


void NonRigidRegistration::operator()( Mesh &mask, const Mesh &target) const
{
    const MatX3f tgtPosVecs = target.features.leftCols(3);
    MatX3f maskPosVecs = mask.features.leftCols(3);

    const K3Tree kdTreeTarget( tgtPosVecs);
    std::shared_ptr<K3Tree> kdTreeMask = std::shared_ptr<K3Tree>( new K3Tree( maskPosVecs));

    const KNNMap kmap( maskPosVecs, *kdTreeMask, _smoothK);
    const SmoothingWeights smw( kmap, mask.flags, _sigmaSmooth);

    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        FlagVec crsFlags; // Return corresponding features on target
        const SparseMat A = _corresponder( *kdTreeMask, mask.flags, kdTreeTarget, target.flags, &crsFlags);
        const FeatMat crs = A * target.features;
        const VecXf iwts = _inlierFinder( mask.features, crs, crsFlags); // Correspondence weights

        // Displacement field from current mask points to correspondence points on target
        DispMat DF = crs.leftCols(3) - mask.features.leftCols(3);

        // Apply visco-elastic steps to the deformation field
        _transformer.update( i, DF, smw, iwts);
        mask.update( DF);

        if ( i < _numUpdateIts - 1)
        {
            // Update for next _corresponder
            maskPosVecs = mask.features.leftCols(3);
            kdTreeMask = std::shared_ptr<K3Tree>( new K3Tree( maskPosVecs));
        }   // end if
    }   // end for
}   // end operator()
