### CMakeLists.txt --- PreCICE library
## Author: Angelos Mantzaflaris, Hugo Verhelst
######################################################################

if(POLICY CMP0076)# CMake 3.13
cmake_policy(SET CMP0076 NEW)
endif()

## PreCICE extension
project(gsPreCICE)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
#aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

find_package(precice REQUIRED CONFIG)

add_library(${PROJECT_NAME} INTERFACE)
# target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)
# target_link_libraries(${PROJECT_NAME} PRIVATE precice::precice)
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY INTERFACE_LINK_LIBRARIES precice::precice)

target_include_directories(${PROJECT_NAME} INTERFACE
 $<INSTALL_INTERFACE:gismo/gsPreCICE>
)

add_dependencies(${PROJECT_NAME} precice)

set(gismo_LINKER ${gismo_LINKER} ${PROJECT_NAME}
   CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "gsPreCICE/*.h")

# Look for CMakeLists.txt in examples
add_subdirectory(examples)