# =============================================================================
#  CADET-semi-analytic - The semi analytic extension of
#			CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2015: Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================
 
# Require a fairly new cmake version
cmake_minimum_required (VERSION 2.6)

# Prohibit in-source build
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message (FATAL_ERROR "In-source build prohibited.")
endif ()

# Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Git)

# Write the current version number to variable and file
if (GIT_FOUND)
  if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE CASEMA_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    # In case of missing tags, default to versions.txt file
    if (NOT "${CASEMA_VERSION}" STREQUAL "")
      message(STATUS "Get version from git")
      
      # Remove first character ("v")
      string(LENGTH "${CASEMA_VERSION}" CASEMA_VERSION_STRLEN)
      math(EXPR CASEMA_VERSION_STRLEN "${CASEMA_VERSION_STRLEN}-1")
      string(SUBSTRING "${CASEMA_VERSION}" 1 ${CASEMA_VERSION_STRLEN}  CASEMA_VERSION)
    endif ()
  endif ()
endif ()

if ("${CASEMA_VERSION}" STREQUAL "")
  message(STATUS "Get version from file")
  file(STRINGS "${CMAKE_SOURCE_DIR}/version.txt" CASEMA_VERSION)
endif ()

message(STATUS "Version is ${CASEMA_VERSION}")

# Get current commit hash from git
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
if (NOT DEFINED GIT_SHA1)
  set(GIT_SHA1 "NO-COMMIT-HASH")
  set(GIT_REFSPEC "NO-REFSPEC")
endif()
message(STATUS "Current git HEAD: ${GIT_REFSPEC} SHA1 ${GIT_SHA1}")

# Name of the current project
project (CASEMA_MASTER CXX)

# Set some variables with info on the CASEMA project
set (PACKAGE_VERSION   "${CASEMA_VERSION}")
set (PACKAGE_BUGREPORT "cadet@fz-juelich.de")
set (PACKAGE_NAME      "CASEMA")
set (PACKAGE_STRING    "CASEMA ${PACKAGE_VERSION}")
set (PACKAGE_TARNAME   "casema-${CASEMA_VERSION}")


# ------------------------------------------------------------------------------------------------------------------------
#   Other configuration options
# ------------------------------------------------------------------------------------------------------------------------
# Option that allows users to build release or debug version
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
#   Try to find modules
# ------------------------------------------------------------------------------------------------------------------------

# Shared components
find_package(HDF5 COMPONENTS CXX C HL)
find_package(GMP)
find_package(MPFR)
find_package(MPC)

# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
#   Options
# ------------------------------------------------------------------------------------------------------------------------

option(BUILD_TESTS "Build test programs" ON)
option(KAHAN_SUMMATION "Use Kahan summation for improved accuracy" ON)
option(USE_FADBAD "Enable FADBAD++ support for algorithmic differentiation" OFF)
set(FADBAD_MAXORDER "200" CACHE STRING "Maximum order of differentiation for FADBAD++")

# Variables for config header
set(CASEMA_KAHAN_SUMMATION ${KAHAN_SUMMATION})
set(CASEMA_USE_FADBAD ${USE_FADBAD})

# Write options to file
configure_file("${CMAKE_SOURCE_DIR}/include/CaSeMaConfig.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/CaSeMaConfig.hpp")

# Add include directories for access to exported header files.
include_directories (${CMAKE_CURRENT_BINARY_DIR})

message (STATUS "Kahan summation: ${CASEMA_KAHAN_SUMMATION}")
message (STATUS "FADBAD++ support: ${CASEMA_USE_FADBAD}")

# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
#   Add selected modules to the build system and add the targets to the list of all targets
# ------------------------------------------------------------------------------------------------------------------------

set (CASEMA_TARGETS)

# Shared components
add_subdirectory (src)

# Tests
if (BUILD_TESTS)
  add_subdirectory (test)
endif()

# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Set properties, definitions, install target etc.
# ------------------------------------------------------------------------------------------------------------------------

# Depending on the used compiler, we enable warnings and OpenMP support: -fopenmp -Wall -Wvla
# Disable line-wrapping in gcc/gxx to support eclipse CDT error parsing: -fmessage-length=0
# Enable C++11 - Standard compiler support: -std=c++0x 
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
message(STATUS "C++ Compiler is ${CXX_COMPILER_NAME} (ID ${CMAKE_CXX_COMPILER_ID}; BIN ${CMAKE_CXX_COMPILER})")
foreach (_TARGET ${CASEMA_TARGETS})

    # Add version info #defines
    set_target_properties (${_TARGET} PROPERTIES COMPILE_DEFINITIONS CASEMA_VERSION="${CASEMA_VERSION}")

    # Set compiler flags for C++11 and OpenMP support
    if ((CXX_COMPILER_NAME STREQUAL "icpc") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR (CXX_COMPILER_NAME STREQUAL "icc"))
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-openmp -w2 -diag-disable remark -std=c++0x -xHOST" LINK_FLAGS "-openmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x -march=native" LINK_FLAGS "")
    elseif (CMAKE_COMPILER_IS_GNUCXX OR CXX_COMPILER_NAME STREQUAL "g++" OR CXX_COMPILER_NAME STREQUAL "c++" OR CXX_COMPILER_NAME STREQUAL "g++.exe")
        if (APPLE)
            set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-fopenmp -Wall -Wvla -std=c++0x -fmessage-length=0" LINK_FLAGS "-fopenmp")
        else()
            set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-fopenmp -Wall -Wvla -std=c++0x -fmessage-length=0 -march=native" LINK_FLAGS "-fopenmp")
        endif()
    elseif (CXX_COMPILER_NAME STREQUAL "cl")
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "/openmp")
    endif ()
endforeach ()

# ------------------------------------------------------------------------------------------------------------------------

