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

#include <Types.h>
#include <cassert>
using rNonRigid::Mesh;
using rNonRigid::FaceMat;
using rNonRigid::MatX6f;
using rNonRigid::MatX3f;
using rNonRigid::Mat4f;
using rNonRigid::MatXf;
using rNonRigid::Vec3f;
using rNonRigid::Vec3i;


namespace {

void updateNorms( MatXf &F, const FaceMat &H)
{
    // Then get the area weighted norms of the faces and add to all associated vertices
    const size_t NF = H.rows();
    for ( size_t i = 0; i < NF; ++i)
    {
        const Vec3i vtxs = H.row(i); // Corresponding vertex rows for this face
        const int a = vtxs[0];
        const int b = vtxs[1];
        const int c = vtxs[2];

        // Get the corresponding updated positions
        const Vec3f vA = F.block<1,3>(a,0);
        const Vec3f vB = F.block<1,3>(b,0);
        const Vec3f vC = F.block<1,3>(c,0);

        // Calculate area weighted triangle norm (magnitude is twice the triangle's area)
        const Vec3f fnrm = (vB - vA).cross( vC - vB);

        // Add back to the associated vertices
        F.block<1,3>(a,3) += fnrm;
        F.block<1,3>(b,3) += fnrm;
        F.block<1,3>(c,3) += fnrm;
    }   // end for

    // Finally renormalise the vertex normals.
    const size_t N = F.rows();
    for ( size_t i = 0; i < N; ++i)
        F.block(i, 3, 1, 3).normalize();  // Normalise in place
}   // end updateNorms

}   // end namespace



void Mesh::update( const MatX3f &D)
{
    // First update vertex positions (first three columns of F) by adding the displacement map
    const size_t N = features.rows();
    assert( D.rows() == N);
    features.leftCols<3>() += D;
    if ( features.cols() >= 6 && topology.rows() > 0)
    {
        features.block(0,3,N,3) = MatX3f::Zero(N, 3);
        updateNorms( features, topology);
    }   // end if
}   // end update


void Mesh::transform( const Mat4f &T)
{
    const size_t N = features.rows();
    MatX4f ft( N, 4);
    ft.leftCols<3>() = features.leftCols<3>();
    ft.rightCols<1>() = VecXf::Ones(N);
    features.leftCols<3>() = (ft * T.transpose()).leftCols<3>();
    if ( features.cols() >= 6)
    {
        const Mat3f R = T.block<3,3>(0,0) / T(3,3); // Rotation submatrix (with possible scaling factor)
        features.block(0,3,N,3) = features.block(0,3,N,3) * R.transpose();
    }   // end if
}   // end transform

