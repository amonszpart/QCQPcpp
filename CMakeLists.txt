project(QCQPcpp)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

######################################################

#remember to set CMAKE_BUILD_TYPE=Release
SET( MOSEK_DIR /home/bontius/workspace/3rdparty/mosek/7/tools/platform/linux64x86)

######################################################

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin )#CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib )#CACHE PATH "Build directory" FORCE)


INCLUDE_DIRECTORIES(
    include
    /usr/include/eigen3/
)

SET( H_LIST
    include/qcqpcpp/optProblem.h
    include/qcqpcpp/mosekOptProblem.h
)
SET( HPP_LIST
    include/qcqpcpp/impl/optProblem.hpp
    include/qcqpcpp/impl/mosekOptProblem.hpp
)

# Mosek
add_library( libmosek64 SHARED IMPORTED )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_RELEASE ${MOSEK_DIR}/bin/libmosek64.so )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_DEBUG ${MOSEK_DIR}/bin/libmosek64.so )
INCLUDE_DIRECTORIES( ${MOSEK_DIR}/h )

## LIB
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${SRC_LIST} ${H_LIST} ${HPP_LIST})
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")
set_target_properties( ${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# LINK
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    libmosek64
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
)
