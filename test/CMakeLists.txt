# =============================================================================
#  CADET-semi-analytic - The semi analytic extension of
#			CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2015-2018: Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project (CASEMA_TEST CXX)

# Set variables defining source and include directories
set (CASEMA_SRC_DIR "${CMAKE_SOURCE_DIR}/test")
set (CASEMA_INC_DIR "${CMAKE_SOURCE_DIR}/include" CACHE INTERNAL "Path to exported libcadet headers")
set (CASEMA_MAIN_SRC_DIR "${CMAKE_SOURCE_DIR}/src")

# Add variable CASEMA_SOURCES with the sources for CASEMA
set (COMMON_SOURCES
     ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp
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

add_executable (testCppAD ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/testCppAD.cpp)
target_link_libraries(testCppAD ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS testCppAD)

# ------------------------------------------------------------------------------------------------------------------------

if (USE_FADBAD)

	add_executable (testFADBAD ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/testFADBAD.cpp)
	target_link_libraries(testFADBAD ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

	list (APPEND CASEMA_TARGETS testFADBAD)

endif()

# ------------------------------------------------------------------------------------------------------------------------

add_executable (extraLab ${CASEMA_SRC_DIR}/extraLab.cpp)
target_link_libraries(extraLab ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS extraLab)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (extraConvLab ${CASEMA_SRC_DIR}/extraConvLab.cpp)
target_link_libraries(extraConvLab ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS extraConvLab)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (extraConsLab ${CASEMA_SRC_DIR}/extraConsLab.cpp)
target_link_libraries(extraConsLab ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS extraConsLab)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (scanmom ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/scanMom.cpp)
target_link_libraries(scanmom ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS scanmom)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (qamarAnaMomTest ${CASEMA_SRC_DIR}/qamarAnaMomTest.cpp)
target_link_libraries(qamarAnaMomTest ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS qamarAnaMomTest)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (qamarMomTest ${CASEMA_SRC_DIR}/qamarMomTest.cpp)
target_link_libraries(qamarMomTest ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS qamarMomTest)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (generalMomTest ${CASEMA_SRC_DIR}/generalMomTest.cpp)
target_link_libraries(generalMomTest ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS generalMomTest)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (inletTest ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/inletTest.cpp)
target_link_libraries(inletTest ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS inletTest)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (fourierCoeff ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/fourierCoeff.cpp)
target_link_libraries(fourierCoeff ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS fourierCoeff)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (fourierError ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/fourierError.cpp)
target_link_libraries(fourierError ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS fourierError)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (scanSol ${COMMON_SOURCES} ${CASEMA_SRC_DIR}/scanSol.cpp)
target_link_libraries(scanSol ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES} ${HDF5_LIBRARIES})

list (APPEND CASEMA_TARGETS scanSol)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (checkCenterMom ${CASEMA_SRC_DIR}/checkCenterMom.cpp)
target_link_libraries(checkCenterMom ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS checkCenterMom)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (testQuadrature ${CASEMA_SRC_DIR}/testQuadrature.cpp)
target_link_libraries(testQuadrature ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS testQuadrature)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (testInverseExpInt ${CASEMA_SRC_DIR}/testInverseExpInt.cpp)
target_link_libraries(testInverseExpInt ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS testInverseExpInt)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (testErrorEstimate ${CASEMA_SRC_DIR}/testErrorEstimate.cpp)
target_link_libraries(testErrorEstimate ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS testErrorEstimate)

# ------------------------------------------------------------------------------------------------------------------------

add_executable (testCSVreader ${CMAKE_SOURCE_DIR}/src/CSVReader.cpp ${CASEMA_SRC_DIR}/testCSVreader.cpp)
target_link_libraries(testCSVreader ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

list (APPEND CASEMA_TARGETS testCSVreader)

# ------------------------------------------------------------------------------------------------------------------------

set (CASEMA_TARGETS "${CASEMA_TARGETS}" PARENT_SCOPE)


# Info message
message (STATUS "Added TEST module")

