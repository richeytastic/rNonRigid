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

#include <FastDeformRegistration.h>
#include <SymmetricCorresponder.h>
#include <SmoothingWeights.h>
#include <InlierFinder.h>
#include <ViscoElasticTransformer.h>
#include <memory>
using rNonRigid::FastDeformRegistration;
using rNonRigid::Mesh;


FastDeformRegistration::FastDeformRegistration( size_t n) : _numUpdateIts(n) {}


void FastDeformRegistration::operator()( Mesh &mask, const Mesh &target) const
{
    //const float RelaxRate = expf( logf( 1.0f/50) / (_numUpdateIts-1));
    const size_t N = mask.features.rows();

    const SymmetricCorresponder corresponder( 7, 0.9f, false);
    std::shared_ptr<K6Tree> kdtTgt = std::shared_ptr<K6Tree>( new K6Tree( target.features));
    std::shared_ptr<K6Tree> kdtMsk = std::shared_ptr<K6Tree>( new K6Tree( mask.features));

    // Create the smoothing map for SmoothKNN neighbours. Values will be between 0 and 1.
    // The smoothing map has a row per vertex with SmoothKNN columns.
    static const size_t SMOOTH_KNN = 80;
    static const float SMOOTH_SIGMA = 3.0f;  // Bigger = greater influence per vertex over neighbours
    const KNNMap kmap( mask.features, *kdtMsk, SMOOTH_KNN);
    const SmoothingWeights swts( kmap, mask.flags, SMOOTH_SIGMA);

    static const float INLIER_KAPPA = 4.0f;
    const InlierFinder inlierFinder( INLIER_KAPPA, true, 10);
   
    for ( size_t i = 0; i < _numUpdateIts; ++i)
    {
        // Get corresponding feature on target for each row (vertex) of mask
        const FeatMat crs = corresponder( *kdtMsk, mask.flags, *kdtTgt, target.flags, nullptr);

        // Update smoothing and inlier weights for later regularisation
        const VecXf iwts = inlierFinder( mask.features, crs, crsFlags);
        const MatXf wts = createWeights( swts, iwts);   // N rows X SMOOTH_KNN columns

        // Vector field to correspondence points weighted by the inliers
        MatX3f VF = crs.leftCols(3) - mask.features.leftCols(3);

        //const size_t numRelax = size_t(roundf( 50.0f * powf( RelaxRate, i)));
        //regularise( VF, wts, swts.indices(), 1); // Regularise the vector field

        // For the initial deformation field, use the *current* mask normals to define the direction
        // of deformation, but scale each vector in proportion to the magnitude of difference to
        // the correspondence point. The delta magnitude is also scaled in inverse proportion to
        // the number of remaining iterations (more iterations giving a more gradual and smoother
        // deformation while allowing for greater sensitivity to local surface detail).
        const float stepScale = 1.0f/(_numUpdateIts-i);
        VecXf M(N); // For each mask vertex, obtain the size of the increment
        for ( size_t j = 0; j < N; ++j)
        {
            const Vec3f &vec = mask.features.block<1,3>(j,3);
            //D.row(j) = stepScale * VF.row(j).dot( vec) * vec;
            M[j] = stepScale * VF.row(j).dot( vec);
        }   // end for

        MatX3f D(N, 3);

        mask.update( D);
        if ( i < _numUpdateIts - 1)
            kdtMsk = std::shared_ptr<K6Tree>( new K6Tree(mask.features));
    }   // end for
}   // end operator()
