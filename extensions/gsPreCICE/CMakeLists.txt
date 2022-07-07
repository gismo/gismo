### CMakeLists.txt --- PreCICE library
## Author: Angelos Mantzaflaris, Hugo Verhelst
#
# The library is header-only, the extra code is downloaded in place of
# the source folder of G+Smo
#
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

message("${PROJECT_NAME}_HEADERS")
message("${${PROJECT_NAME}_HEADERS}")

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

# add example files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})
