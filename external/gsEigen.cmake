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

find_package(Eigen3 3.3 NO_MODULE)

if(NOT Eigen3_FOUND)

    message(STATUS "Eigen3 not found, fetching ...")

    include(FetchContent)

    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0)

    # FetchContent_MakeAvailable(eigen EXCLUDE_FROM_ALL)

    FetchContent_GetProperties(eigen)

    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()

    add_library(Eigen3 INTERFACE)

    set(EIGEN3_INCLUDE_DIRS ${eigen_SOURCE_DIR})
endif()

include_directories(${EIGEN3_INCLUDE_DIRS})

message(STATUS "Eigen3 include dirs: ${EIGEN3_INCLUDE_DIRS}")

## #############################################################################
## Code ends here
