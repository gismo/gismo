
if(USE_VCXPROJ)

message(STATUS "Importing MSVC project THBSplines.")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.rc")                              

include_external_msproject(THBSplines
                           ${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.vcxproj
                           GismoSupport)
set_target_properties(THBSplines PROPERTIES
    FOLDER "Rhino THB plug-ins"  )

#---- Native cmake
else()

cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

if(POLICY CMP0043)# CMake 3.1
cmake_policy(SET CMP0043 OLD)
endif()

project(THBSplines CXX RC)

set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

include_directories(${GISMO_INCLUDE_DIRS})

set(${PROJECT_NAME}_HEADERS
    Resource.h
    stdafx.h
    targetver.h
    THBSplinesApp.h
    THBSplinesPlugIn.h
    )
    
set(${PROJECT_NAME}_SOURCES
    cmdActiveLevel.cpp
    cmdFitThbSurface.cpp
    cmdRefineNurbsAtSpan.cpp
    cmdTHBSplines.cpp
    cmdRefineNurbs.cpp
    cmdToNurbs.cpp
    cmdToPatches.cpp
    stdafx.cpp
    THBSplinesApp.cpp
    THBSplinesPlugIn.cpp
    )
  
#set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1) #cmake>=3.4, maybe problematic

include_directories(${RHINO_SDK}/inc)

add_definitions(-D_UNICODE -DUNICODE -D_USRDLL -DWIN64 -UWIN32 -DDO_EXPORTS) # -DRHINO_V6_READY

add_library(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_HEADERS}
  ${${PROJECT_NAME}_SOURCES}
  ${PROJECT_NAME}.rc
  )

set_target_properties (${PROJECT_NAME} PROPERTIES
        COMPILE_DEFINITIONS _AFXDLL;_BIND_TO_CURRENT_CRT_VERSION;_BIND_TO_CURRENT_MFC_VERSION
        #LINK_FLAGS "/ENTRY:\"wWinMainCRTStartup\""
)

# set "link library dependencies" to true
#set_target_properties(${PROJECT_NAME} PROPERTIES STATIC_LIBRARY_FLAGS "\" LinkLibraryDependencies=\"true")

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(${PROJECT_NAME}  
#${RHINO_SDK}/lib/Release/rhino.lib
#${RHINO_SDK}/lib/Release/rdk.lib
#${RHINO_SDK}/lib/Release/opennurbs.lib
GismoSupport
)

target_link_libraries(${PROJECT_NAME} gismo ${gismo_LINKER} )

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plRhino)

install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION include/gismo/plugins
  LIBRARY DESTINATION include/gismo/plugins
  ARCHIVE DESTINATION include/gismo/plugins)

# This target is a plugin, so we hook it to the grouping target
add_dependencies(plugins ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME} PROPERTIES
    FOLDER "Rhino THB plug-ins"  )

endif() #native cmake

