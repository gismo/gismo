######################################################################
### CMakeLists.txt --- gsParasolid
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

## Parasolid extension
project(gsParasolidExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)
#message ("*** ${PROJECT_NAME} sources: ${${PROJECT_NAME}_HEADERS}, ${${PROJECT_NAME}_SOURCES}")

# Look for parasolid library
find_package(Parasolid REQUIRED)
include_directories(${PARASOLID_INCLUDE_DIR})

# Set default schema path  
add_definitions(-DGISMO_P_SCHEMA="${PARASOLID_INCLUDE_DIR}/schema")

# Compile statically the extension
add_library(${PROJECT_NAME} OBJECT
${${PROJECT_NAME}_SOURCES}
${${PROJECT_NAME}_HEADERS}
${${PROJECT_NAME}_HPPFILES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "${PROJECT_NAME} extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${PARASOLID_LIBRARY}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

if(GISMO_BUILD_LIB)
install(FILES   ${CMAKE_CURRENT_SOURCE_DIR}/gsReadParasolid.h
		${CMAKE_CURRENT_SOURCE_DIR}/gsWriteParasolid.h
		${CMAKE_CURRENT_SOURCE_DIR}/gsPKSession.h
        DESTINATION include/gismo/gsParasolid/)
endif()

#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
