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
project(CasemaTest CXX C)

add_executable(testRunner testRunner.cpp
	BesselZeros.cpp
	DurbinsMethod.cpp
	DinisExpansion.cpp
	ExpInt.cpp
	EstimateTests.cpp
	GRM2D.cpp
	SystemTests.cpp
	$<TARGET_OBJECTS:casema_common>)

target_link_libraries(testRunner PRIVATE CASEMA::CompileOptions GMP::GMP MPFR::MPFR MPC::MPC Eigen3::Eigen)
if (OpenMP_FOUND)
	target_link_libraries(testRunner PRIVATE OpenMP::OpenMP_CXX)
endif()
target_include_directories(testRunner PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR})
target_include_directories(testRunner PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/Catch ${CMAKE_SOURCE_DIR}/ThirdParty/mpreal ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml)

# Info message
message(STATUS "Added tests")

