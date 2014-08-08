project(QCQPcpp)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

include_directories( include )

SET( H_LIST
    include/qcqpcpp/optProblem.h
    include/qcqpcpp/mosekOptProblem.h
)
SET( HPP_LIST
    include/qcqpcpp/impl/optProblem.hpp
    include/qcqpcpp/impl/mosekOptProblem.hpp
)

# Mosek
SET( MOSEK_DIR /home/bontius/workspace/3rdparty/mosek/7/tools/platform/linux64x86 )
add_library( libmosek64 SHARED IMPORTED )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_RELEASE ${MOSEK_DIR}/bin/libmosek64.so )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_DEBUG ${MOSEK_DIR}/bin/libmosek64.so )
INCLUDE_DIRECTORIES( ${MOSEK_DIR}/h )

# EXECUTABLE
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${SRC_LIST} ${H_LIST} ${HPP_LIST})
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

# LINK
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    libmosek64
)
