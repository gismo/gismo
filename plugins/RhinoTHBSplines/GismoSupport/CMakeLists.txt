
if(true)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.rc") # sets version value

include_external_msproject (GismoSupport
                           ${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.vcxproj )

set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

#---- Native cmake
else()

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

if(POLICY CMP0043)# CMake 3.1
cmake_policy(SET CMP0043 OLD)
endif()

project(GismoSupport CXX RC)

set(${PROJECT_NAME}_HEADERS
    ActiveLevelCommandImplementation.h
    CommandImplementation.h
    FileExportImplementation.h
    FileImportImplementation.h
    GismoExportImplementation.h
    GismoImportImplementation.h
    gsStreamData.h
    RefineNurbsAtSpanCommandImplementation.h
    RhinoGetThbSurfaceObject.h
    TestCommandImplementation.h
    ThbControlNetDrawer.h
    ThbFitCommandImplementation.h
    ThbGripsState.h
    ThbSurfaceDrawer.h
    ThbSurfaceGrip.h
    ThbSurfaceGrips.h
    ThbSurfaceObject.h
    ThbSurfaceUserData.h
    ThbSurfaceUtils.h
    ThbSurfaceWatcher.h
    RefineNurbsCommandImplementation.h
    Resource.h
    stdafx.h
    targetver.h
    ToNurbsImplementation.h
    ToPatchesImplementation.h
)

set(${PROJECT_NAME}_SOURCES
    ActiveLevelCommandImplementation.cpp
    GismoExportImplementation.cpp
    GismoImportImplementation.cpp
    gsStreamData.cpp
    RefineNurbsAtSpanCommandImplementation.cpp
    RhinoGetThbSurfaceObject.cpp
    TestCommandImplementation.cpp
    ThbControlNetDrawer.cpp
    ThbFitCommandImplementation.cpp
    ThbGripsState.cpp
    ThbSurfaceDrawer.cpp
    ThbSurfaceGrip.cpp
    ThbSurfaceGrips.cpp
    ThbSurfaceObject.cpp
    ThbSurfaceUserData.cpp
    ThbSurfaceUtils.cpp
    ThbSurfaceWatcher.cpp
    RefineNurbsCommandImplementation.cpp
    stdafx.cpp
    ToNurbsImplementation.cpp
    ToPatchesImplementation.cpp
    GismoSupport.rc
    )

add_definitions(-D_AFXDLL)
set(CMAKE_MFC_FLAG 2)

add_library(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_HEADERS}
  ${${PROJECT_NAME}_SOURCES}
  GismoSupport.rc
  )

target_link_libraries(${PROJECT_NAME} gismo_static ${gismo_LINKER} )

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plRhino)

install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION include/gismo/plugins
  LIBRARY DESTINATION include/gismo/plugins
  ARCHIVE DESTINATION include/gismo/plugins)

# This target is a plugin, so we hook it to the grouping target
add_dependencies(plugins ${PROJECT_NAME})

endif() #native cmake


set_target_properties(${PROJECT_NAME} PROPERTIES
    FOLDER "Rhino THB plug-ins"  )
