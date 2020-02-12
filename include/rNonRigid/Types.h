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

#ifndef RNONRIGID_TYPES_H
#define RNONRIGID_TYPES_H

#include "rNonRigid_Export.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#ifdef _WIN32
// Disable warning about DLL linkage to Eigen library not being exported (it's header only).
#pragma warning( disable : 4251)
#endif

namespace rNonRigid {

static const int NFEATURES = 6;

using FeatMat = Eigen::Matrix<float, Eigen::Dynamic, NFEATURES>; // N x 6 features matrix
using FaceMat = Eigen::Matrix<int, Eigen::Dynamic, 3>;      // M x 3 integers giving vertices indices per face
using DispMat = Eigen::Matrix<float, Eigen::Dynamic, 3>;    // Displacement field of vector deltas (xyz)
using MatX3f = DispMat;

using FeatVec = Eigen::Matrix<float, NFEATURES, 1>;         // A single feature in NFEATURES space
using FlagVec = Eigen::VectorXf;                    // Dynamic length vector of floats with values either 1 or 0

using Vec3f = Eigen::Vector3f;                      // Point in 3 space
using Mat3f = Eigen::Matrix3f;                      // 3x3 matrix (for cross variance)
using Vec4f = Eigen::Vector4f;                      // Point in 4 space
using Mat4f = Eigen::Matrix4f;                      // 4x4 matrix (for transforms)
using VecXf = Eigen::VectorXf;                      // Dynamic length vector of floats with values in [0,1]
using MatXf = Eigen::MatrixXf;                      // Dynamic size matrix of floats
using MatXi = Eigen::MatrixXi;                      // Dynamic size matrix of integers

using SparseMat = Eigen::SparseMatrix<float>;

struct rNonRigid_EXPORT Mesh
{
    FeatMat features;   // Features (vertices and normals) per row
    FaceMat topology;   // Face connectivity as row indices into features
    FlagVec flags;      // Features (vertices) to use as 1's or 0's (flags.size() == features.rows())
    void update( const DispMat&);   // Calls rNonRigid::updateFeatures
};  // end Mesh

/**
 * Add DispMat to the first three columns (position) of FeatMat, then
 * update vertex normals in last three columns using the topology.
 */
rNonRigid_EXPORT void updateFeatures( FeatMat&, const FaceMat&, const DispMat&);

}   // end namespace

#endif
