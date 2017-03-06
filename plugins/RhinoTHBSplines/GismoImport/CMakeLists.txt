configure_file("${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GismoImport.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/GismoImport.rc")
                              

include_external_msproject(GismoImport
                           ${CMAKE_CURRENT_SOURCE_DIR}/GismoImport.vcxproj
						   GismoSupport)
						   
set_target_properties(GismoImport PROPERTIES
    FOLDER "Rhino THB plug-ins"  )