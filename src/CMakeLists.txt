# =============================================================================
#  CADET-semi-analytic - The semi analytic extension of
#			CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2015-2017: Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project (CASEMA CXX)

# Set variables defining source and include directories
set (CASEMA_SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set (CASEMA_INC_DIR "${CMAKE_SOURCE_DIR}/include" CACHE INTERNAL "Path to headers")

# Git integration
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/VersionInfo.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp" @ONLY)

# Add variable CASEMA_SOURCES with the sources for CASEMA
set (COMMON_SOURCES
     ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp
     ${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp
    )

# Add include directories for access to exported header files.
include_directories (${CASEMA_INC_DIR} ${CASEMA_SRC_DIR})

# Add lib headers
include_directories (${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR} ${MPC_INCLUDE_DIR})

# Add ThirdParty headers
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/mpreal")
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/pugixml")
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include")
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/quadpackpp/include")
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/cppad")
if (USE_FADBAD)
	include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/FADBAD++")
endif()
# ------------------------------------------------------------------------------------------------------------------------

add_executable (chrom ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/chrom.cpp)
target_link_libraries(chrom ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS chrom)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (moments ${COMMON_SOURCES} ${CMAKE_SOURCE_DIR}/src/CSVReader.cpp ${CASEMA_SRC_DIR}/moments.cpp)
target_link_libraries(moments ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS moments)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (singlePointSolution ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/singlePointSolution.cpp)
target_link_libraries(singlePointSolution ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS singlePointSolution)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (singlePointMethodComp ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/singlePointMethodComp.cpp)
target_link_libraries(singlePointMethodComp ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS singlePointMethodComp)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (momentQuadrature ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/momentQuadrature.cpp)
target_link_libraries(momentQuadrature ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS momentQuadrature)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (chromError ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/chromError.cpp ${CASEMA_SRC_DIR}/CSVReader.cpp)
target_link_libraries(chromError ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS chromError)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (errorCalc ${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp ${CASEMA_SRC_DIR}/errorCalc.cpp ${CASEMA_SRC_DIR}/CSVReader.cpp)
target_link_libraries(errorCalc ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS errorCalc)

# ------------------------------------------------------------------------------------------------------------------------

set (CASEMA_TARGETS "${CASEMA_TARGETS}" PARENT_SCOPE)


# Installing
install (CODE "MESSAGE(\"\nInstall CASEMA\n\")")
install (TARGETS ${CASEMA_TARGETS} DESTINATION bin)


# Info message
message (STATUS "Added CASEMA module")

