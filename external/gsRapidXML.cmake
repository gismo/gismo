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

find_package(RapidXML 1.13)

if(NOT RapidXML_FOUND)

    message(STATUS "RapidXML not found, fetching ...")

    include(FetchContent)

    FetchContent_Declare(
        rapidxml
        URL https://sourceforge.net/projects/rapidxml/files/rapidxml/rapidxml%201.13/rapidxml-1.13.zip/download
        URL_HASH SHA1=f5fd4fbc5ad7e96045313697811d65ea8089a950
    )

    FetchContent_Populate(rapidxml)

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${rapidxml_SOURCE_DIR} ${CMAKE_BINARY_DIR}/rapidxml
        RESULT_VARIABLE result
    )

    set(RapidXML_INCLUDE_DIR ${CMAKE_BINARY_DIR})
endif()

include_directories(${RapidXML_INCLUDE_DIR})

message(STATUS "RapidXML include dir: ${RapidXML_INCLUDE_DIR}")

## #############################################################################
## Code ends here
