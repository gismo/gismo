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

find_package(opennurbs QUIET)

if(NOT opennurbs_FOUND)

    message(STATUS "Opennurbs not found, fetching ...")

    include(FetchContent)
    include(GNUInstallDirs)

    FetchContent_Declare(
        opennurbs
        GIT_REPOSITORY https://gitlab.inria.fr/gismo/opennurbs.git
        GIT_TAG master
    )

    FetchContent_MakeAvailable(opennurbs)

    set(OPENNURBS_INCLUDE_DIR ${opennurbs_SOURCE_DIR})
endif()

include_directories(${OPENNURBS_INCLUDE_DIR})

## #############################################################################
## Code ends here
