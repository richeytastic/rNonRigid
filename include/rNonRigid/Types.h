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

using FaceMat = Eigen::Matrix<int, Eigen::Dynamic, 3>;  // M x 3 integers of vertex indices per face
using MatX6f = Eigen::Matrix<float, Eigen::Dynamic, 6>; // N x 6 features matrix
using MatX4f = Eigen::Matrix<float, Eigen::Dynamic, 4>; // For homogeneous matrix multiplication
using MatX3f = Eigen::Matrix<float, Eigen::Dynamic, 3>; // Displacement field of vector deltas (xyz)

using FeatVec = Eigen::Matrix<float, 6, 1>;         // A single feature

using Vec3i = Eigen::Vector3i;
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
    Mesh() {}
    explicit Mesh( size_t n) : features(n,6) {}

    MatX6f features;   // Features (vertices and normals) per row
    FaceMat topology;  // Face connectivity as row indices into features

    void refreshNormals();

    void update( const MatX3f&);    // Calls rNonRigid::updateFeatures

    // Return a set of features by adding given displacement field to this mesh's features.
    MatX6f makeFeatures( const MatX3f&) const;
};  // end Mesh

/**
 * Add MatX3f to the first three columns (position) of MatX6f, then
 * update vertex normals in last three columns using the topology.
 */
rNonRigid_EXPORT void updateFeatures( MatX6f&, const FaceMat&, const MatX3f&);

/**
 * Transform matrix M by T where the first three columns of M should
 * be the X,Y,Z coordinates of the row positions. The remaining
 * columns can be anything invariant to affine transform.
 */
rNonRigid_EXPORT void transform( MatXf &M, const Mat4f& T);

}   // end namespace

#endif
