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

# Name of the current project
project(CasemaApps CXX C)

# Git integration
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/VersionInfo.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp" ESCAPE_QUOTES @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/CompileTimeConfig.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/CompileTimeConfig.hpp" @ONLY)

# Compiler features
write_compiler_detection_header(
	FILE "${CMAKE_BINARY_DIR}/casemaCompilerInfo.hpp"
	PREFIX CASEMA
	COMPILERS GNU Clang AppleClang MSVC Intel
	FEATURES cxx_noexcept cxx_user_literals cxx_constexpr cxx_variadic_templates
	PROLOG "\
// =============================================================================\n\
//  CADET-semi-analytic - The semi-analytic extension of CADET\n\
//  \n\
//  Copyright © 2015-2020: Samuel Leweke¹²\n\
//  \n\
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.\n\
//    ² University of Cologne, Cologne, Germany.\n\
//  \n\
//  All rights reserved. This program and the accompanying materials\n\
//  are made available under the terms of the GNU Public License v3.0 (or, at\n\
//  your option, any later version) which accompanies this distribution, and\n\
//  is available at http://www.gnu.org/licenses/gpl.html\n\
// =============================================================================\n"
	EPILOG "\
#if CASEMA_COMPILER_CXX_CONSTEXPR\n\
	#define CASEMA_CONST_OR_CONSTEXPR constexpr\n\
#else\n\
	#define CASEMA_CONST_OR_CONSTEXPR const\n\
#endif\n\
#if CASEMA_COMPILER_CXX_USER_LITERALS && CASEMA_COMPILER_CXX_CONSTEXPR\n\
	#define CASEMA_COMPILETIME_HASH 1\n\
#else\n\
	#define CASEMA_COMPILETIME_HASH 0\n\
#endif\n"
)

# Check for some compiler bugs
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state(RESET)

# Enable C++11
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	if (WIN32)
		set(CMAKE_REQUIRED_FLAGS "/Qstd=c++11")
	else()
		set(CMAKE_REQUIRED_FLAGS "-std=c++11")
	endif()
elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
	set(CMAKE_REQUIRED_FLAGS "-std=c++11")
elseif ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_REQUIRED_FLAGS "-std=c++11")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	set(CMAKE_REQUIRED_FLAGS "")
endif ()

set(CMAKE_REQUIRED_INCLUDES "${CMAKE_BINARY_DIR}/")
CHECK_CXX_SOURCE_COMPILES("#include \"casemaCompilerInfo.hpp\" \n class Test { public: Test() { }\n Test& operator=(Test&& cpy) CASEMA_NOEXCEPT = default; \n private: \n int a; }; int main(int argc, char** argv) { return 0; }" COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE)

cmake_pop_check_state()

# Collect common stuff
add_library(casema_common OBJECT "${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp" 
	GraphAlgos.cpp
	BesselZeros.cpp
	ModelBuilder.cpp
	TimeReader.cpp
	DurbinsMethod.cpp
	LaplaceError.cpp
	model/UnitOperation.cpp
	model/ConvectionDispersionLikeModel.cpp
	model/ColumnLikeModel.cpp
	model/LumpedRateModelWithoutPores.cpp
	model/LumpedRateModelWithPores.cpp
	model/GeneralRateModel.cpp
	model/GeneralRateModel2D.cpp
	model/StirredTankModel.cpp
	model/InletModel.cpp
	model/OutletModel.cpp
	model/ModelSystem.cpp)

target_link_libraries(casema_common PUBLIC CASEMA::CompileOptions GMP::GMP MPFR::MPFR MPC::MPC Eigen3::Eigen)

target_include_directories(casema_common PUBLIC 
	"${CMAKE_BINARY_DIR}"
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_SOURCE_DIR}/include/"
	"${CMAKE_SOURCE_DIR}/ThirdParty/json"
	"${CMAKE_SOURCE_DIR}/ThirdParty/pugixml"
	"${CMAKE_SOURCE_DIR}/ThirdParty/mpreal")



add_library(casema_io OBJECT ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp
	io/FileIO.cpp)

target_link_libraries(casema_io PUBLIC CASEMA::CompileOptions)

target_include_directories(casema_io PUBLIC 
	"${CMAKE_BINARY_DIR}"
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_SOURCE_DIR}/include/"
	"${CMAKE_SOURCE_DIR}/ThirdParty/json"
	"${CMAKE_SOURCE_DIR}/ThirdParty/pugixml")




add_executable(chromatogram $<TARGET_OBJECTS:casema_common> $<TARGET_OBJECTS:casema_io> ProgressBar.cpp chromatogram.cpp)
target_link_libraries(chromatogram PUBLIC CASEMA::CompileOptions PRIVATE GMP::GMP MPFR::MPFR MPC::MPC Eigen3::Eigen)
target_include_directories(chromatogram PUBLIC 
	"${CMAKE_BINARY_DIR}"
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_SOURCE_DIR}/include/"
	"${CMAKE_SOURCE_DIR}/ThirdParty/json"
	"${CMAKE_SOURCE_DIR}/ThirdParty/pugixml"
	"${CMAKE_SOURCE_DIR}/ThirdParty/mpreal"
	"${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include")

if (OpenMP_FOUND)
	target_link_libraries(casema_common PUBLIC OpenMP::OpenMP_CXX)
	target_link_libraries(chromatogram PUBLIC OpenMP::OpenMP_CXX)
endif()
if (HDF5_FOUND)
	target_link_libraries(casema_io PUBLIC HDF5::HDF5)
	target_link_libraries(chromatogram PUBLIC HDF5::HDF5)
endif()

# ---------------------------------------------------
#   Set installation related stuff
# ---------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CHROMATOGRAM\n\")")
install(TARGETS chromatogram RUNTIME)

## Info message
message(STATUS "Added CHROMATOGRAM module")
