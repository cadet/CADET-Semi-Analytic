# =============================================================================
#  CADET-semi-analytic - The semi-analytic extension of CADET
#  
#  Copyright © 2015-2020: Samuel Leweke¹²
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#    ² University of Cologne, Cologne, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================
 
# Require a fairly new cmake version
cmake_minimum_required(VERSION 3.14)

# Prohibit in-source build
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
	message(FATAL_ERROR "In-source build prohibited.")
endif()

# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Git)

# Write the current version number to variable
if (GIT_FOUND)
	if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
		execute_process(COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 HEAD
		                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		                OUTPUT_VARIABLE CASEMA_VERSION
		                OUTPUT_STRIP_TRAILING_WHITESPACE)

		if (NOT "${CASEMA_VERSION}" STREQUAL "")
			message(STATUS "Get version from git")

			# Remove first character ("v")
			string(LENGTH "${CASEMA_VERSION}" CASEMA_VERSION_STRLEN)
			math(EXPR CASEMA_VERSION_STRLEN "${CASEMA_VERSION_STRLEN}-1")
			string(SUBSTRING "${CASEMA_VERSION}" 1 ${CASEMA_VERSION_STRLEN}  CASEMA_VERSION)
		endif()
	endif()
endif()

# In case of missing tags, default to versions.txt file
if ("${CASEMA_VERSION}" STREQUAL "")
	message(STATUS "Get version from file")
	file(STRINGS "${CMAKE_SOURCE_DIR}/version.txt" CASEMA_VERSION)
endif()

message(STATUS "CASEMA version: ${CASEMA_VERSION}")

# Get current commit hash from git
if (GIT_FOUND)
	include(GetGitRevisionDescription)
	get_git_head_revision(GIT_REFSPEC GIT_SHA1)
endif()
if (NOT DEFINED GIT_SHA1)
	set(GIT_SHA1 "NO-COMMIT-HASH")
	set(GIT_REFSPEC "NO-REFSPEC")
endif()
message(STATUS "Current git HEAD: ${GIT_REFSPEC} SHA1 ${GIT_SHA1}")

# Name of the current project
project(CaSeMaFramework
	VERSION ${CASEMA_VERSION}
	DESCRIPTION "Liquid column chromatography reference solutions"
	HOMEPAGE_URL "https://github.com/modsim/CADET-semi-analytic"
	LANGUAGES CXX C)

# Always use '-fPIC'/'-fPIE' option
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Hide symbols by default
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# Enable folders for IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ---------------------------------------------------
#   Other configuration options
# ---------------------------------------------------

include(FeatureSummary)

option(ENABLE_PLATFORM_TIMER "Use a platform-dependent timer" OFF)
add_feature_info(ENABLE_PLATFORM_TIMER ENABLE_PLATFORM_TIMER "Use a platform-dependent timer")

option(ENABLE_BESSEL_TRUNCATION "Truncate output of Bessel functions near 0.0 to 0.0" OFF)
add_feature_info(ENABLE_BESSEL_TRUNCATION ENABLE_BESSEL_TRUNCATION "Truncate output of Bessel functions near 0.0 to 0.0")

option(ENABLE_TESTS "Build CASEMA tests" ON)
add_feature_info(ENABLE_TESTS ENABLE_TESTS "Build CASEMA tests")


option(ENABLE_IPO "Enable interprocedural optimization if compiler supports it" ON)
add_feature_info(ENABLE_IPO ENABLE_IPO "Enable interprocedural optimization if compiler supports it")

option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath" ON)
add_feature_info(CMAKE_INSTALL_RPATH_USE_LINK_PATH CMAKE_INSTALL_RPATH_USE_LINK_PATH "Add paths to linker search and installed rpath")

# Hande RPATH on OSX when not installing to a system directory, see
# https://groups.google.com/d/msg/fenics-dev/KSCrob4M_1M/zsJwdN-SCAAJ
# and https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
if (APPLE)
	# The RPATH to be used when installing, but only if it's not a system directory
	set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
	if ("${isSystemDir}" STREQUAL "-1")
		set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	endif()
endif()

# Option that allows users to build release or debug version
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (default)")
endif()

# ---------------------------------------------------
#   Check build environment
# ---------------------------------------------------
include(WriteCompilerDetectionHeader)

find_package(GMP)
find_package(MPFR)
find_package(MPC)
find_package(OpenMP COMPONENTS CXX)
find_package(Eigen3)

#if (ON)
	set(HDF5_USE_STATIC_LIBRARIES ${ENABLE_STATIC_LINK_DEPS})
	find_package(HDF5 COMPONENTS C)
	set_package_properties(HDF5 PROPERTIES
		DESCRIPTION "Hierarchical Data Format 5 (HDF5)"
		URL "https://www.hdfgroup.org/HDF5"
		TYPE RECOMMENDED
		PURPOSE "File IO"
	)

	if (HDF5_FOUND)
		add_library(HDF5::HDF5 UNKNOWN IMPORTED)

		list(LENGTH HDF5_C_LIBRARIES HDF5_C_LEN)
		if (HDF5_C_LEN GREATER 1)
			list(GET HDF5_C_LIBRARIES 0 HDF5_MAIN_LIBRARY)
			set(HDF5_SUPPORT_LIBRARIES ${HDF5_C_LIBRARIES})
			list(REMOVE_AT HDF5_SUPPORT_LIBRARIES 0)
		else()
			set(HDF5_MAIN_LIBRARY ${HDF5_C_LIBRARIES}) 
			set(HDF5_SUPPORT_LIBRARIES)
		endif()

		set_target_properties(HDF5::HDF5 PROPERTIES
			IMPORTED_LOCATION ${HDF5_MAIN_LIBRARY}
			INTERFACE_INCLUDE_DIRECTORIES "${HDF5_C_INCLUDE_DIRS}"
#			INTERFACE_COMPILE_DEFINITIONS ${HDF5_C_DEFINITIONS}
		)
		if (HDF5_SUPPORT_LIBRARIES)
			target_link_libraries(HDF5::HDF5 INTERFACE ${HDF5_SUPPORT_LIBRARIES})
		endif()

		unset(HDF5_SUPPORT_LIBRARIES)
		unset(HDF5_MAIN_LIBRARY)
		unset(HDF5_C_LEN)

		# Make sure HDF5_LIBRARY_DIRS is defined
		if ((NOT DEFINED HDF5_LIBRARY_DIRS) OR (NOT HDF5_LIBRARY_DIRS) OR ("${HDF5_LIBRARY_DIRS}" STREQUAL ""))
			list(GET HDF5_LIBRARIES 0 HDF5_LIB_TEMP)
			get_filename_component(HDF5_LIBRARY_DIRS ${HDF5_LIB_TEMP} DIRECTORY)
			unset(HDF5_LIB_TEMP)
		endif()

		# Check if we need additional libraries for linking (i.e., zlib, szip)
		include(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
		include(${CMAKE_ROOT}/Modules/CMakePushCheckState.cmake)

		cmake_push_check_state(RESET)

		# Set libs and includes
		set(CMAKE_REQUIRED_LIBRARIES ${HDF5_LIBRARIES})
		set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS})

		CHECK_CXX_SOURCE_COMPILES("#include <hdf5.h>\nint main(int argc, char** argv){\n H5Zfilter_avail(H5Z_FILTER_SZIP);\nH5Zfilter_avail(H5Z_FILTER_DEFLATE);\nreturn 0;\n}\n" HDF5_DONT_NEED_ZLIBS)

		# Reset libs and includes
		cmake_pop_check_state()

		# Find szip and zlib libs if we need them
		if (NOT HDF5_DONT_NEED_ZLIBS)

			# Prefer static libs if enabled
			set(_HDF5_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
			if(ENABLE_STATIC_LINK_DEPS)
				if(WIN32)
					set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
				else()
					set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
				endif()
			endif()

			find_library(HDF5_SZLIB NAMES libszip szip PATHS ${HDF5_LIBRARY_DIRS})
			find_library(HDF5_ZLIB NAMES libzlib zlib z libz PATHS ${HDF5_LIBRARY_DIRS})

			if (HDF5_SZLIB)
				list(APPEND HDF5_LIBRARIES ${HDF5_SZLIB})
				add_library(HDF5::SZLIB UNKNOWN IMPORTED)
				set_target_properties(HDF5::SZLIB PROPERTIES IMPORTED_LOCATION ${HDF5_SZLIB})
				target_link_libraries(HDF5::HDF5 INTERFACE HDF5::SZLIB)
			endif()
			if (HDF5_ZLIB)
				list(APPEND HDF5_LIBRARIES ${HDF5_ZLIB})
				add_library(HDF5::ZLIB UNKNOWN IMPORTED)
				set_target_properties(HDF5::ZLIB PROPERTIES IMPORTED_LOCATION ${HDF5_ZLIB})
				target_link_libraries(HDF5::HDF5 INTERFACE HDF5::ZLIB)
			endif()
			unset(HDF5_SZLIB)
			unset(HDF5_ZLIB)

			set(CMAKE_FIND_LIBRARY_SUFFIXES ${_HDF5_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
			unset(_HDF5_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
		endif()
	endif()
#endif()


set(ENABLE_HDF5 "${HDF5_FOUND}")
set(IPO_AVAILABLE OFF)
if (ENABLE_IPO)
	include(CheckIPOSupported)
	check_ipo_supported(RESULT IPO_RESULT OUTPUT IPO_OUT LANGUAGES CXX)
	if (IPO_RESULT)
		set(IPO_AVAILABLE ON)
		set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
	else()
		message(WARNING "IPO is not supported: ${IPO_OUT}")
	endif()
	unset(IPO_RESULT)
	unset(IPO_OUT)
endif()


# ---------------------------------------------------
#   Add selected modules to the build system and add the targets to the list of all targets
# ---------------------------------------------------

add_library(CASEMA::CompileOptions INTERFACE IMPORTED)
target_compile_features(CASEMA::CompileOptions INTERFACE cxx_std_14)
set(CMAKE_CXX_EXTENSIONS OFF)

if (WIN32)
	target_compile_definitions(CASEMA::CompileOptions INTERFACE NOMINMAX)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_definitions(CASEMA::CompileOptions INTERFACE CASEMA_LOGLEVEL_MIN=Trace DEBUG _DEBUG)
else()
	target_compile_definitions(CASEMA::CompileOptions INTERFACE CASEMA_LOGLEVEL_MIN=Warning NDEBUG)
endif()

if (ENABLE_PLATFORM_TIMER)
	target_compile_definitions(CASEMA::CompileOptions INTERFACE CASEMA_USE_PLATFORM_TIMER)
	if ((NOT APPLE) AND (NOT WIN32))
		target_link_libraries(CASEMA::CompileOptions INTERFACE rt)
	endif()
endif()

# CASEMA library
add_subdirectory(src)

if (ENABLE_TESTS)
	add_subdirectory(test)
endif()


# ---------------------------------------------------
#   Set properties, definitions, install target etc.
# ---------------------------------------------------

# Packaging support
set(CPACK_PACKAGE_VENDOR "CASEMA")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
	/.git
	/\\\\.DS_Store
)

message("")
message("--------------------------- Feature Summary ---------------------------")
         
feature_summary(WHAT ALL)

# Summary
message("")
message("------------------------------- Summary -------------------------------")
message("C++ compiler name: ${CMAKE_CXX_COMPILER_ID} at ${CMAKE_CXX_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
message("Source dir: ${CMAKE_SOURCE_DIR}")
message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install dir: ${CMAKE_INSTALL_PREFIX}")
message("C Flags: ${CMAKE_C_FLAGS}")
message("C++ Flags: ${CMAKE_CXX_FLAGS}")
message("IPO enabled: ${IPO_AVAILABLE}")
message("------------------------------- Modules -------------------------------")
message("Tests: ${ENABLE_TESTS}")
message("------------------------------- Options -------------------------------")
message("Platform-dependent timer: ${ENABLE_PLATFORM_TIMER}")
message("Truncation of Bessel functions: ${ENABLE_BESSEL_TRUNCATION}")
message("----------------------------- Dependencies ----------------------------")

message("Found Eigen: ${EIGEN3_FOUND}")
if (EIGEN3_FOUND)
	message("  Version ${EIGEN3_VERSION}")
	message("  Include ${EIGEN3_INCLUDE_DIR}")
endif()

message("Found GMP: ${GMP_FOUND}")
if (GMP_FOUND)
	message("  Version ${GMP_VERSION}")
	message("  Include ${GMP_INCLUDE_DIR}")
	message("  Libs ${GMP_LIBRARIES}")
endif()

message("Found MPFR: ${MPFR_FOUND}")
if (MPFR_FOUND)
	message("  Version ${MPFR_VERSION}")
	message("  Include ${MPFR_INCLUDE_DIR}")
	message("  Libs ${MPFR_LIBRARIES}")
endif()

message("Found MPC: ${MPC_FOUND}")
if (MPC_FOUND)
	message("  Version ${MPC_VERSION}")
	message("  Include ${MPC_INCLUDE_DIR}")
	message("  Libs ${MPC_LIBRARIES}")
endif()

message("Found OpenMP: ${OpenMP_FOUND}")
if (OpenMP_FOUND)
	message("  Version ${OpenMP_VERSION}")
	message("  Include ${OpenMP_CXX_INCLUDE_DIRS}")
	message("  Compile flags ${OpenMP_CXX_FLAGS}")
	message("  Libs ${OpenMP_CXX_LIBRARIES}")
endif()

message("Found HDF5: ${HDF5_FOUND}")
if (HDF5_FOUND)
	message("  Version ${HDF5_VERSION}")
	message("  Includes ${HDF5_INCLUDE_DIRS}")
	message("  Libs ${HDF5_LIBRARIES}")
	message("  Defs ${HDF5_C_DEFINITIONS}")
endif()

message("-----------------------------------------------------------------------")
message("")
