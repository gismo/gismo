 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

## #################################################################
## Collect Headers
## #################################################################

project(gsParaview)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} HFILES)
aux_cxx_source_directory(${CMAKE_CURRENT_SOURCE_DIR} CFILES)
#message ("*** Paraview sources: ${HFILES}, ${CFILES}")

# Locate Qt
FIND_PACKAGE( Qt4 REQUIRED )
QT4_WRAP_CPP(MOC_SRCS ${HFILES})
QT4_ADD_RESOURCES(RCS_SRCS GismoGeometryPanel.qrc)

# Locate ParaView build and import CMake configuration, 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
 
############################
###### Add paraview objects
############################

ADD_PARAVIEW_ACTION_GROUP(IFACES IFACE_SRCS
                            CLASS_NAME SourceToolbarActions
                            GROUP_NAME "ToolBar/SourceToolbar"
)

ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS CLASS_NAME GismoGeometryPanel
                            XML_NAME GismoSource XML_GROUP sources)


############################
###### Build Paraview plugin
############################

 ADD_PARAVIEW_PLUGIN(
   GismoParaview               #<--Name for the plugin
  "0.1"                        #<--Version string
  SERVER_MANAGER_XML gsParaview.xml
   SOURCES ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS} ${CFILES} 
   SERVER_MANAGER_SOURCES GismoGeometrySource.cxx
   GUI_INTERFACES ${IFACES}
   GUI_RESOURCE_FILES pqReader.xml
)


set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/plParaview)
add_custom_target( gsParaview DEPENDS GismoParaview)

#target_link_libraries( GismoParaview gismo)



#
# http://www.paraview.org/Wiki/ParaView/Users_Guide/Loading_Plugins
#
#
