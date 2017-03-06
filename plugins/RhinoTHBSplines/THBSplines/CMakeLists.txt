configure_file("${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/PluginVersion.h")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.rc")                              

include_external_msproject(THBSplines
                           ${CMAKE_CURRENT_SOURCE_DIR}/THBSplines.vcxproj
						   GismoSupport)
set_target_properties(THBSplines PROPERTIES
    FOLDER "Rhino THB plug-ins"  )	