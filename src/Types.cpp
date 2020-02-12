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

#include <Types.h>
#include <cassert>
using rNonRigid::Mesh;
using rNonRigid::FeatMat;
using rNonRigid::FaceMat;
using rNonRigid::DispMat;


void rNonRigid::updateFeatures( FeatMat &F, const FaceMat &H, const DispMat &D)
{
    // First update vertex positions (first three columns of F) by adding the displacement map
    const size_t N = F.rows();
    assert( int(N) == D.rows());
    for ( size_t i = 0; i < N; ++i)
    {
        F.block<1,3>(i,0) += D.row(i);
        F.block<1,3>(i,3) = Vec3f::Zero();  // Zero out the norms (recalculating)
    }   // end for

    // Then get the area weighted norms of the faces and add to all associated vertices
    const size_t NF = H.rows();
    for ( size_t i = 0; i < NF; ++i)
    {
        const Eigen::Vector3i vtxs = H.row(i); // Corresponding vertex rows for this face

        // Get the corresponding updated positions
        const Vec3f vA = F.block<1,3>( vtxs[0], 0);
        const Vec3f vB = F.block<1,3>( vtxs[1], 0);
        const Vec3f vC = F.block<1,3>( vtxs[2], 0);

        // Calculate the area weighted triangle norm
        const Vec3f fnrm = (vB - vA).cross( vC - vB);  // DON'T NORMALISE (magnitude is twice triangle's area)

        // Add back to the associated vertices
        F.block<1,3>(vtxs[0], 3) += fnrm;
        F.block<1,3>(vtxs[1], 3) += fnrm;
        F.block<1,3>(vtxs[2], 3) += fnrm;
    }   // end for

    // Finally renormalise the vertex normals.
    for ( size_t i = 0; i < N; ++i)
        F.block<1,3>(i, 3).normalize();  // Normalise in place
}   // end updateFeatures


void Mesh::update( const DispMat &D) { updateFeatures( features, topology, D);}
