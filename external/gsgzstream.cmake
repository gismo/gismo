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

find_package(gzstream QUIET)

if(NOT gzstream_FOUND)

    message(STATUS "gzstream not found, fetching ...")

    include(FetchContent)

    FetchContent_Declare(
        gzstream
        GIT_REPOSITORY https://gitlab.inria.fr/gismo/gzstream.git
        GIT_TAG main
    )

    FetchContent_MakeAvailable(gzstream)

    message(STATUS "gzstream source dir: ${gzstream_SOURCE_DIR}")

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${gzstream_SOURCE_DIR} ${CMAKE_BINARY_DIR}/gzstream
        RESULT_VARIABLE result
    )

endif()

## #############################################################################
## Code ends here
