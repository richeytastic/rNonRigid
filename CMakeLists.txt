cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

project( rNonRigid)

set(WITH_EIGEN TRUE)
set(WITH_NANOFLANN TRUE)
include( "cmake/FindLibs.cmake")

set( INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set( INCLUDE_F "${INCLUDE_DIR}/${PROJECT_NAME}")
set( SRC_DIR "${PROJECT_SOURCE_DIR}/src")
include_directories( ${INCLUDE_DIR})
include_directories( ${INCLUDE_F})

set( INCLUDE_FILES
    "${INCLUDE_F}.h"
    "${INCLUDE_F}/Types.h"
    "${INCLUDE_F}/InlierFinder.h"
    #"${INCLUDE_F}/FastDeformRegistration.h"
    "${INCLUDE_F}/K6Tree.h"
    "${INCLUDE_F}/K3Tree.h"
    "${INCLUDE_F}/KNNCorresponder.h"
    "${INCLUDE_F}/KNNMap.h"
    "${INCLUDE_F}/NonRigidRegistration.h"
    "${INCLUDE_F}/NonSymmetricCorresponder.h"
    "${INCLUDE_F}/RigidRegistration.h"
    "${INCLUDE_F}/RigidTransformer.h"
    "${INCLUDE_F}/SmoothingWeights.h"
    "${INCLUDE_F}/SymmetricCorresponder.h"
    "${INCLUDE_F}/ViscoElasticTransformer.h"
    )

set( SRC_FILES
    ${SRC_DIR}/InlierFinder
    #${SRC_DIR}/FastDeformRegistration
    ${SRC_DIR}/K6Tree
    ${SRC_DIR}/K3Tree
    ${SRC_DIR}/KNNCorresponder
    ${SRC_DIR}/KNNMap
    ${SRC_DIR}/NonRigidRegistration
    ${SRC_DIR}/NonSymmetricCorresponder
    ${SRC_DIR}/RigidRegistration
    ${SRC_DIR}/RigidTransformer
    ${SRC_DIR}/SmoothingWeights
    ${SRC_DIR}/SymmetricCorresponder
    ${SRC_DIR}/Types
    ${SRC_DIR}/ViscoElasticTransformer
    )

add_library( ${PROJECT_NAME} ${SRC_FILES} ${INCLUDE_FILES})
include( "cmake/LinkLibs.cmake")
