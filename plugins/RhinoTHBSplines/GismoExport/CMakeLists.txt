configure_file("${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GismoExport.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/GismoExport.rc")
                              

include_external_msproject(GismoExport
                           ${CMAKE_CURRENT_SOURCE_DIR}/GismoExport.vcxproj
						   GismoSupport)
						   
set_target_properties(GismoExport PROPERTIES
    FOLDER "Rhino THB plug-ins"  )