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

void updateNorms( MatX6f &F, const FaceMat &H)
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
        const Vec3f vA = F.block<1,3>( a, 0);
        const Vec3f vB = F.block<1,3>( b, 0);
        const Vec3f vC = F.block<1,3>( c, 0);

        // Calculate area weighted triangle norm (magnitude is twice the triangle's area)
        const Vec3f fnrm = (vB - vA).cross( vC - vB);

        // Add back to the associated vertices
        F.block<1,3>(a, 3) += fnrm;
        F.block<1,3>(b, 3) += fnrm;
        F.block<1,3>(c, 3) += fnrm;
    }   // end for

    // Finally renormalise the vertex normals.
    const size_t N = F.rows();
    for ( size_t i = 0; i < N; ++i)
        F.block<1,3>(i, 3).normalize();  // Normalise in place
}   // end updateNorms

}   // end namespace


void rNonRigid::updateFeatures( MatX6f &F, const FaceMat &H, const MatX3f &D)
{
    // First update vertex positions (first three columns of F) by adding the displacement map
    assert( D.rows() == F.rows());
    F.leftCols(3) += D;
    F.rightCols(3) = MatX3f::Zero(F.rows(), 3);
    updateNorms( F, H);
}   // end updateFeatures


void Mesh::update( const MatX3f &D) { updateFeatures( features, topology, D);}


void Mesh::refreshNormals()
{
    features.rightCols(3) = MatX3f::Zero(features.rows(), 3);
    updateNorms( features, topology);
}   // end refreshNormals


MatX6f Mesh::makeFeatures( const MatX3f &D) const
{
    assert( features.rows() == D.rows());
    MatX6f nf( features.rows(), 6);
    nf.leftCols(3) = features.leftCols(3) + D;
    nf.rightCols(3) = Mat3f::Zero( features.rows(), 3);
    updateNorms( nf, topology);
    return nf;
}   // end makeFeatures


void rNonRigid::transform( MatXf &F, const Mat4f &T)
{
    MatX4f ft( F.rows(), 4);
    ft.leftCols(3) = F.leftCols(3);
    ft.rightCols(1) = VecXf::Ones(ft.rows());
    ft = ft * T.transpose();
    F.leftCols(3) = ft.leftCols(3);
}   // end transform

