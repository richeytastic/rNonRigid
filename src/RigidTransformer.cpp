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

#include <RigidTransformer.h>
#include <iostream>
#include <cassert>
using RNonRigid::RigidTransformer;
using RNonRigid::FeatMat;
using RNonRigid::Vec3f;
using RNonRigid::Vec4f;
using RNonRigid::VecXf;
using RNonRigid::Mat3f;
using RNonRigid::Mat4f;
using RNonRigid::MatXf;

namespace {

Vec3f findWeightedMean( const MatXf& ps, const VecXf& wts, float wsum)
{
    const size_t N = wts.size();
    assert( N == size_t(ps.cols()));
    Vec3f wm = Vec3f::Zero();
    for ( size_t i = 0; i < N; ++i)
        wm += wts[i] * ps.col(i);
    return wm / wsum;
}   // end findWeightedMean


Mat3f computeCrossVariance( const MatXf& fps, const Vec3f& fwm,
                            const MatXf& cps, const Vec3f& cwm,
                            const VecXf& wts, float wsum)
{
    const size_t N = wts.size();
    assert( N == size_t(fps.cols()));
    assert( N == size_t(cps.cols()));
    Mat3f cmat = Mat3f::Zero();
    for ( size_t i = 0; i < N; ++i)
        cmat += wts[i] * fps.col(i) * cps.col(i).transpose();
    return cmat / wsum - fwm * cwm.transpose();
}   // end computeCrossVariance


Mat4f computeQ( const Mat3f& cv)
{
    const Mat3f asymm = cv - cv.transpose();   // Anti-symmetric matrix
    const Vec3f delta( asymm(1,2), asymm(2,0), asymm(0,1));    // Use cyclic elements to construct delta
    const float trace = cv.trace();
    Mat4f Q;
    Q(0,0) = trace;
    Q.block<1,3>(0,1) = delta;
    Q.block<3,1>(1,0) = delta;
    Q.block<3,3>(1,1) = cv + cv.transpose() - trace * Mat3f::Identity();
    return Q;
}   // end computeQ


Mat3f computeRotationMatrix( const Mat4f& Q)
{
    Eigen::SelfAdjointEigenSolver<Mat4f> decomposer(Q);
    assert( decomposer.info() == Eigen::Success);
    size_t idxMaxVal = 0;
    float maxEigenVal = decomposer.eigenvalues()[0];
    for ( size_t i = 1; i < 4; ++i)
    {
        if ( decomposer.eigenvalues()[i] > maxEigenVal)
        {
            maxEigenVal = decomposer.eigenvalues()[i];
            idxMaxVal = i;
        }   // end if
    }   // end for

    // Rotation quaternion as eigenvector of Q with largest eigenvalue
    const Vec4f rq = decomposer.eigenvectors().col( idxMaxVal);

    Mat3f R;
    const float rq0s = rq(0) * rq(0);
    const float rq1s = rq(1) * rq(1);
    const float rq2s = rq(2) * rq(2);
    const float rq3s = rq(3) * rq(3);

    R(0,0) = rq0s + rq1s - rq2s - rq3s;
    R(1,1) = rq0s + rq2s - rq1s - rq3s;
    R(2,2) = rq0s + rq3s - rq1s - rq2s;

    R(0,1) = 2.0f * (rq(1) * rq(2) - rq(0) * rq(3));
    R(1,0) = 2.0f * (rq(1) * rq(2) + rq(0) * rq(3));
    R(0,2) = 2.0f * (rq(1) * rq(3) + rq(0) * rq(2));
    R(2,0) = 2.0f * (rq(1) * rq(3) - rq(0) * rq(2));
    R(1,2) = 2.0f * (rq(2) * rq(3) - rq(0) * rq(1));
    R(2,1) = 2.0f * (rq(2) * rq(3) + rq(0) * rq(1));

    return R;
}   // end computeRotationMatrix


float estimateScaleFactor( const MatXf& fps, const Vec3f& fwm,
                           const MatXf& cps, const Vec3f& cwm,
                           const VecXf& wts, const Mat3f& R)
{
    const size_t N = wts.size();
    assert( N == size_t(fps.cols()));
    assert( N == size_t(cps.cols()));
    float num = 0.0f;
    float den = 0.0f;
    for ( size_t i = 0; i < N; ++i)
    {
        const Vec3f nfp = R * (fps.col(i) - fwm);   // Centre and rotate the floating position
        const Vec3f ncp = cps.col(i) - cwm;         // Centre the corresponding position
        num += wts[i] * ncp.dot(nfp);
        den += wts[i] * nfp.dot(nfp);
    }   // end for
    assert( den != 0.0);
    return num/den;
}   // end estimateScaleFactor

}   // end namespace


RigidTransformer::RigidTransformer( bool us) : _useScaling(us) {}


Mat4f RigidTransformer::operator()( const FeatMat &flt, const FeatMat &crs, const VecXf &wts) const
{
    assert( flt.rows() == crs.rows());
    assert( flt.rows() == wts.size());
    assert( NFEATURES == crs.cols());
    assert( NFEATURES == flt.cols());

    // Set the position triplets to be column vectors (ignore the remaining columns)
    MatXf fps = flt.leftCols(3).transpose();
    MatXf cps = crs.leftCols(3).transpose();

    // Find the weighted mean of the positions for each set
    const float wsum = wts.sum();
    const Vec3f fwm = findWeightedMean( fps, wts, wsum);
    const Vec3f cwm = findWeightedMean( cps, wts, wsum);

    // Get the rotation, scaling, and translation components
    const Mat3f CV = computeCrossVariance( fps, fwm, cps, cwm, wts, wsum);
    const Mat4f Q = computeQ( CV);

    const Mat3f R = computeRotationMatrix( Q);
    const float sf = _useScaling ? estimateScaleFactor( fps, fwm, cps, cwm, wts, R) : 1.0f; // Scale factor
    const Vec3f t = cwm - sf * R * fwm; // Translation between centroids

    // Construct final transformation matrix for return
    Mat4f T = Mat4f::Zero();
    T.block<3,1>(0,3) = t;       // Translation
    T.block<3,3>(0,0) = sf * R;  // Rotation
    T(3,3) = sf;                 // Scale factor in bottom right
    return T;
}   // end operator()


// Apply T to the features which should be 6 elements long (columns) in each row with the
// first three elements giving a point's position, and the last three elements giving its normal.
FeatMat RNonRigid::transformFeatures( const FeatMat& F, const Mat4f& T)
{
    assert( F.cols() == NFEATURES);
    const size_t N = F.rows();
    FeatMat W( N, NFEATURES);

    const float sf = T(3,3);    // Scale factor in bottom right corner
    Mat4f R = Mat4f::Identity();
    R.block<3,3>(0,0) = T.block<3,3>(0,0) / sf;  // Rotation submatrix

    Vec4f ipos = Vec4f::Ones();
    Vec4f inrm = Vec4f::Ones();
    for ( size_t i = 0; i < N; ++i)
    {
        ipos.head<3>() = F.block<1,3>(i,0);
        const Vec4f opos = T * ipos;  // Transform position
        W.block<1,3>(i,0) = opos.head<3>();

        inrm.head<3>() = F.block<1,3>(i,3);
        const Vec4f onrm = R * inrm;  // Apply rotation matrix to normal and remove scale factor
        W.block<1,3>(i,3) = onrm.head<3>();
    }   // end for

    return W;
}   // end transformFeatures
