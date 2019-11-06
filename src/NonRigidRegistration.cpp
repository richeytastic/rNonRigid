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

#include <NonRigidRegistration.h>
#include <SmoothingWeights.h>
#include <iostream>
#include <cassert>
#include <ctime>
using RNonRigid::NonRigidRegistration;
using RNonRigid::Mat4f;
using RNonRigid::Vec3f;
using RNonRigid::FeatMat;
using RNonRigid::DispMat;
using RNonRigid::FaceMat;
using RNonRigid::FlagVec;


void NonRigidRegistration::updateFeatures( FeatMat &uF, const FaceMat& H, const DispMat &D)
{
    // First update vertex positions (first three columns of uF) by adding the displacement map
    const size_t N = uF.rows();
    assert( int(N) == D.rows());
    for ( size_t i = 0; i < N; ++i)
    {
        uF.block<1,3>(i,0) += D.row(i);
        uF.block<1,3>(i,3) = Vec3f::Zero();  // Zero out the norms (recalculating)
    }   // end for

    // Then get the area weighted norms of the faces and add to all associated vertices
    const size_t NF = H.rows();
    for ( size_t i = 0; i < NF; ++i)
    {
        const Eigen::Vector3i vtxs = H.row(i); // Corresponding vertex rows for this face

        // Get the corresponding updated positions
        const Vec3f vA = uF.block<1,3>( vtxs[0], 0);
        const Vec3f vB = uF.block<1,3>( vtxs[1], 0);
        const Vec3f vC = uF.block<1,3>( vtxs[2], 0);

        // Calculate the area weighted triangle norm
        const Vec3f fnrm = (vB - vA).cross( vC - vB);  // DON'T NORMALISE (magnitude is twice triangle's area)

        // Add back to the associated vertices
        uF.block<1,3>(vtxs[0], 3) += fnrm;
        uF.block<1,3>(vtxs[1], 3) += fnrm;
        uF.block<1,3>(vtxs[2], 3) += fnrm;
    }   // end for

    // Finally renormalise the vertex normals.
    for ( size_t i = 0; i < N; ++i)
        uF.block<1,3>(i, 3).normalize();  // Normalise in place
}   // end updateFeatures


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


void NonRigidRegistration::operator()( FeatMat &uF, const FaceMat &H, const FlagVec &fF, const FeatMat &tgt, const FlagVec &tF) const
{
    //Profiling::Profiler profiler( "RNonRigid::NonRigidRegistration::operator()", false);

    //profiler.startSplit("Initialise");
    const KNNMap kmap( uF, uF, _smoothK);
    const SmoothingWeights smw( kmap, fF, _sigmaSmooth);
    const KDTree tkdt( tgt);
    //profiler.endSplit("Initialise");

    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        //profiler.startSplit("Correspondences");
        FlagVec cFlags;
        const FeatMat crs = _corresponder( uF, fF, tkdt, tF, cFlags);   // Returns corresponding features on target
        //profiler.endSplit("Correspondences");

        //profiler.startSplit("Inlier Finding");
        const VecXf wts = _inlierFinder( uF, crs, cFlags);              // Returns correspondence weightings
        //profiler.endSplit("Inlier Finding");

        //profiler.startSplit("Transforming");
        const DispMat CF = crs.leftCols(3) - uF.leftCols(3);            // Raw force field from current floating to corresponding points
        const DispMat D = _transformer.update( i, CF, smw, wts, fF);    // Calculate the displacement for this iteration
        updateFeatures( uF, H, D);
        //profiler.endSplit("Transforming");
    }   // end for
}   // end operator()
