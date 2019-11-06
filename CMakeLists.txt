cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project( rNonRigid)

set(WITH_EIGEN TRUE)
set(WITH_NANOFLANN TRUE)
include( "cmake/FindLibs.cmake")

set( INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set( SRC_DIR "${PROJECT_SOURCE_DIR}/src")
include_directories( ${INCLUDE_DIR})

set( INCLUDE_FILES
    "${INCLUDE_DIR}/rNonRigid.h"
    "${INCLUDE_DIR}/Types.h"
    "${INCLUDE_DIR}/InlierFinder.h"
    "${INCLUDE_DIR}/FeaturesCalculate.h"
    "${INCLUDE_DIR}/KDTree.h"
    "${INCLUDE_DIR}/KNNCorresponder.h"
    "${INCLUDE_DIR}/KNNMap.h"
    "${INCLUDE_DIR}/NonRigidRegistration.h"
    "${INCLUDE_DIR}/RigidRegistration.h"
    "${INCLUDE_DIR}/RigidTransformer.h"
    "${INCLUDE_DIR}/SmoothingWeights.h"
    "${INCLUDE_DIR}/SymmetricCorresponder.h"
    "${INCLUDE_DIR}/ViscoElasticTransformer.h"
    )

set( SRC_FILES
    ${SRC_DIR}/InlierFinder
    ${SRC_DIR}/FeaturesCalculate
    ${SRC_DIR}/KDTree
    ${SRC_DIR}/KNNCorresponder
    ${SRC_DIR}/KNNMap
    ${SRC_DIR}/NonRigidRegistration
    ${SRC_DIR}/RigidRegistration
    ${SRC_DIR}/RigidTransformer
    ${SRC_DIR}/SmoothingWeights
    ${SRC_DIR}/SymmetricCorresponder
    ${SRC_DIR}/ViscoElasticTransformer
    )

add_library( ${PROJECT_NAME} ${SRC_FILES} ${INCLUDE_FILES})
include( "cmake/LinkLibs.cmake")
