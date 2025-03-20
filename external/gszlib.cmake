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

set(ZLIB_USE_STATIC_LIBS ON)

find_package(ZLIB)

if(NOT ZLIB_FOUND)
    message(STATUS "ZLIB not found, fetching ...")

    include(FetchContent)

    FetchContent_Declare(
        zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1)

    FetchContent_MakeAvailable(zlib)

    set(ZLIB_INCLUDE_DIRS ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
else()
    add_library(zlibstatic ALIAS ZLIB::ZLIB)
endif()

include_directories(${ZLIB_INCLUDE_DIRS})

## #############################################################################
## Code ends here
