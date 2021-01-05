# - Try to find libIRIS_C
# Once done this will define
#  IRIS_C_FOUND - System has libIRIS_C
#  IRIS_C_INCLUDE_DIRS - The libIRIS_C include directories
#  IRIS_C_LIBRARIES - The libraries needed to use libIRIS_C

FIND_PATH(WITH_IRIS_C_PREFIX
    NAMES include/iris_c.h
    HINTS /usr/gapps/iris/current/$ENV{SYS_TYPE}/
)

FIND_LIBRARY(IRIS_C_LIBRARIES
    NAMES iris_c
    HINTS ${WITH_IRIS_C_PREFIX}/lib
)
IF(NOT IRIS_C_LIBRARIES)
   FIND_LIBRARY(IRIS_C_LIBRARIES
        NAMES iris_c
        HINTS ${WITH_IRIS_C_PREFIX}/lib
    )
ENDIF(NOT IRIS_C_LIBRARIES) 

FIND_PATH(IRIS_C_INCLUDE_DIRS
    NAMES iris_c.h
    HINTS ${WITH_IRIS_C_PREFIX}/include
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(IRIS DEFAULT_MSG
    IRIS_C_LIBRARIES
    IRIS_C_INCLUDE_DIRS
)

# Hide these vars from ccmake GUI
MARK_AS_ADVANCED(
	IRIS_C_LIBRARIES
	IRIS_C_INCLUDE_DIRS
)