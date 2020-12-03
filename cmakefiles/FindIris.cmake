# - Try to find libIRIS
# Once done this will define
#  IRIS_FOUND - System has libIRIS
#  IRIS_INCLUDE_DIRS - The libIRIS include directories
#  IRIS_LIBRARIES - The libraries needed to use libIRIS

FIND_PATH(WITH_IRIS_PREFIX
    NAMES include/iris.hpp
    HINTS /usr/gapps/iris/current/$ENV{SYS_TYPE}/
)

FIND_LIBRARY(IRIS_LIBRARIES
    NAMES iris
    HINTS ${WITH_IRIS_PREFIX}/lib
)
IF(NOT IRIS_LIBRARIES)
   FIND_LIBRARY(IRIS_LIBRARIES
        NAMES iris
        HINTS ${WITH_IRIS_PREFIX}/lib
    )
ENDIF(NOT IRIS_LIBRARIES) 

FIND_PATH(IRIS_INCLUDE_DIRS
    NAMES iris.hpp
    HINTS ${WITH_IRIS_PREFIX}/include
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(IRIS DEFAULT_MSG
    IRIS_LIBRARIES
    IRIS_INCLUDE_DIRS
)

# Hide these vars from ccmake GUI
MARK_AS_ADVANCED(
	IRIS_LIBRARIES
	IRIS_INCLUDE_DIRS
)