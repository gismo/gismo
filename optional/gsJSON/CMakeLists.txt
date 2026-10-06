### CMakeLists.txt --- JSON library
## Author: H.M. Verhelst
#
# The library is header-only, the extra code is downloaded in place of
# the source folder of G+Smo
#
######################################################################

## JSON extension
project(gsJSONExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)

# Apply same configuration as G+Smo
include(gsConfig)

include(gsFetch)
gismo_fetch_directory(JSON
  URL https://github.com/nlohmann/json/releases/download/v3.12.0/include.zip
  DESTINATION external
)

# Turn off JSON's test
set(JSON_BuildTests OFF CACHE INTERNAL "")

set(JSON_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/JSON/include/ CACHE INTERNAL "")
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${JSON_INCLUDE_DIR}
    CACHE INTERNAL "Gismo include directories" FORCE)

include_directories(${GISMO_INCLUDE_DIRS})

add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)
add_dependencies(${PROJECT_NAME} JSON)

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  #START Export all symbols from this extension
  CXX_VISIBILITY_PRESET default
  C_VISIBILITY_PRESET default
  VISIBILITY_INLINES_HIDDEN 0
  #END Export all symbols from this extension
  FOLDER "G+Smo extensions"
  )

# For gsXmlUtils.h
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/JSON/include/nlohmann
  DESTINATION include/gismo
  PATTERN "*/.svn" EXCLUDE
  )

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "gsJSON/*.h")

# add example files
#include_directories(${SPECTRA_INCLUDE_DIR})
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})


# add unittests
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
