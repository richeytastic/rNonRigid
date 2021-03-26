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

#include <RigidTransformer.h>
#include <cassert>
#include <cfloat>
using rNonRigid::RigidTransformer;
using rNonRigid::Vec3f;
using rNonRigid::Vec4f;
using rNonRigid::VecXf;
using rNonRigid::Mat3f;
using rNonRigid::Mat4f;
using rNonRigid::MatXf;

namespace {

// Returns cross variance matrix (vertices stored per row for fps/cps)
Mat3f computeCV( const MatXf& fps, const Vec3f& fwm,
                  const MatXf& cps, const Vec3f& cwm,
                  const VecXf& wts, float wsum)
{
    const MatXf wfps = fps.array().colwise() * wts.array();
    return ((wfps.transpose() * cps) / wsum) - fwm * cwm.transpose();
}   // end computeCV


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


// For vectors stored as rows in fps/cps
float estimateScaleFactor( const MatXf& fps, const Vec3f& fwm,
                            const MatXf& cps, const Vec3f& cwm,
                            const VecXf& wts, const Mat3f& R)
{
    const MatXf cfp = (fps.rowwise() - fwm.transpose()) * R.transpose(); // Centre and rotate and weight the floating positions
    const MatXf ccp = cps.rowwise() - cwm.transpose(); // Centre the corresponding positions
    const float den = cfp.cwiseProduct( cfp).rowwise().sum().dot(wts) + FLT_MIN;
    const float num = ccp.cwiseProduct( cfp).rowwise().sum().dot(wts) + FLT_MIN;
    return num/den;
}   // end estimateScaleFactor

}   // end namespace


RigidTransformer::RigidTransformer( bool us) : _useScaling(us) {}


Mat4f RigidTransformer::operator()( const MatX3f &flt, const MatX3f &crs, const VecXf &wts) const
{
    assert( flt.rows() == crs.rows());
    assert( flt.rows() == wts.size());

    // Find the weighted mean of the positions for each set
    const float wsum = wts.sum();
    const Vec3f fwm = (wts.transpose() * flt) / wsum;
    const Vec3f cwm = (wts.transpose() * crs) / wsum;

    // Get the rotation, scaling, and translation components
    const Mat3f R = computeRotationMatrix( computeQ( computeCV( flt, fwm, crs, cwm, wts, wsum)));
    const float sf = _useScaling ? estimateScaleFactor( flt, fwm, crs, cwm, wts, R) : 1.0f; // Scale factor
    const Vec3f t = cwm - sf * R * fwm; // Translation between centroids

    // Construct final transformation matrix for return
    Mat4f T = Mat4f::Zero();
    T.block<3,1>(0,3) = t;       // Translation
    T.block<3,3>(0,0) = sf * R;  // Rotation
    T(3,3) = sf;                 // Scale factor in bottom right
    return T;
}   // end operator()
