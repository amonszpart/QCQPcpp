project(QCQPcpp)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

######################################################

#remember to set CMAKE_BUILD_TYPE=Release
SET( MOSEK_DIR "/home/bontius/workspace/3rdparty/mosek/7/tools/platform/linux64x86" )
SET( BONMIN_DIR "/home/bontius/workspace/3rdparty/CoinBonmin-stable/build/")

######################################################

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin )#CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib )#CACHE PATH "Build directory" FORCE)


##_____________________________________________________________________________
# Mosek
add_library( libmosek64 SHARED IMPORTED )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_RELEASE ${MOSEK_DIR}/bin/libmosek64.so )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_DEBUG ${MOSEK_DIR}/bin/libmosek64.so )
INCLUDE_DIRECTORIES( ${MOSEK_DIR}/h )

##_____________________________________________________________________________
# Bonmin
SET( BONMIN_INCLUDE_DIRS
    ${BONMIN_DIR}/include
)
SET( BONMIN_LIBRARIES
    ${BONMIN_DIR}/lib/libipopt.so
    ${BONMIN_DIR}lib/libbonmin.so
    ${BONMIN_DIR}lib/libCbc.so
    #${BONMIN_DIR}lib/libCbcSolver.so
    ${BONMIN_DIR}lib/libCgl.so
    ${BONMIN_DIR}lib/libClp.so
    ${BONMIN_DIR}lib/libcoinhsl.so
    ${BONMIN_DIR}lib/libCoinUtils.so
    ${BONMIN_DIR}lib/libOsi.so
    #${BONMIN_DIR}lib/libOsiCbc.so
    ${BONMIN_DIR}lib/libOsiClp.so
    #${BONMIN_DIR}lib/libOsiCommonTests.so
    lapack
    dl
)

INCLUDE_DIRECTORIES(
    include
    /usr/include/eigen3/
    ${BONMIN_INCLUDE_DIRS}
)

SET( H_LIST
    include/qcqpcpp/optProblem.h
    include/qcqpcpp/mosekOptProblem.h
    include/qcqpcpp/bonminOptProblem.h
)
SET( HPP_LIST
    include/qcqpcpp/impl/optProblem.hpp
    include/qcqpcpp/impl/mosekOptProblem.hpp
)

##_____________________________________________________________________________
## LIB
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${SRC_LIST} ${H_LIST} ${HPP_LIST})
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")
set_target_properties( ${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# LINK
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    libmosek64
    ${BONMIN_LIBRARIES}
)

## TEST
SET( TEST_SRC_LIST
    src/test/test.cpp
)

SET( TEST_PROJECT_NAME "test_${PROJECT_NAME}" )
add_executable( ${TEST_PROJECT_NAME} ${H_LIST} ${HPP_LIST} ${TEST_SRC_LIST} )
set_target_properties(${TEST_PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# LINK
TARGET_LINK_LIBRARIES( ${TEST_PROJECT_NAME}
    libmosek64
    ${BONMIN_LIBRARIES}
)
