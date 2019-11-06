# =============================================================================
# The rNonRigid CMake configuration file.
#
#           ** File generated automatically, DO NOT MODIFY! ***

# To use from an external project, in your project's CMakeLists.txt add:
#   FIND_PACKAGE( rNonRigid REQUIRED)
#   INCLUDE_DIRECTORIES( rNonRigid ${rNonRigid_INCLUDE_DIRS})
#   LINK_DIRECTORIES( ${rNonRigid_LIBRARY_DIR})
#   TARGET_LINK_LIBRARIES( MY_TARGET_NAME ${rNonRigid_LIBRARIES})
#
# This module defines the following variables:
#   - rNonRigid_FOUND         : True if rNonRigid is found.
#   - rNonRigid_ROOT_DIR      : The root directory where rNonRigid is installed.
#   - rNonRigid_INCLUDE_DIRS  : The rNonRigid include directories.
#   - rNonRigid_LIBRARY_DIR   : The rNonRigid library directory.
#   - rNonRigid_LIBRARIES     : The rNonRigid imported libraries to link to.
#
# =============================================================================

get_filename_component( rNonRigid_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component( rNonRigid_ROOT_DIR  "${rNonRigid_CMAKE_DIR}"           PATH)

set( rNonRigid_INCLUDE_DIRS "${rNonRigid_ROOT_DIR}/../include" CACHE PATH "The rNonRigid include directories.")
set( rNonRigid_LIBRARY_DIR  "${rNonRigid_ROOT_DIR}"            CACHE PATH "The rNonRigid library directory.")

include( "${CMAKE_CURRENT_LIST_DIR}/Macros.cmake")
get_library_suffix( _lsuff)
set( _hints rNonRigid${_lsuff} librNonRigid${_lsuff})
find_library( rNonRigid_LIBRARIES NAMES ${_hints} PATHS "${rNonRigid_LIBRARY_DIR}/static" "${rNonRigid_LIBRARY_DIR}")
set( rNonRigid_LIBRARIES     ${rNonRigid_LIBRARIES}         CACHE FILEPATH "The rNonRigid imported libraries to link to.")

# handle QUIETLY and REQUIRED args and set rNonRigid_FOUND to TRUE if all listed variables are TRUE
include( "${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake")
find_package_handle_standard_args( rNonRigid rNonRigid_FOUND rNonRigid_LIBRARIES rNonRigid_INCLUDE_DIRS)
