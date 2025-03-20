## Version: $Id:  $
##
##

## Commentary:
##
##

## Changelog:
##
##

##
## Code starts here
## #############################################################################

find_package(TCLAP 1.4 QUIET)

if(NOT TCLAP_FOUND)

    message(STATUS "TCLAP not found, fetching ...")

    include(FetchContent)

    FetchContent_Declare(
        tclap
        GIT_REPOSITORY https://github.com/mirror/tclap.git
        GIT_TAG 1.4
    )

    FetchContent_Populate(tclap)

    set(TCLAP_INCLUDE_DIRS ${tclap_SOURCE_DIR}/include)
endif()

include_directories(${TCLAP_INCLUDE_DIRS})

message(STATUS "TCLAP include dirs: ${TCLAP_INCLUDE_DIRS}")

## #############################################################################
## Code ends here
