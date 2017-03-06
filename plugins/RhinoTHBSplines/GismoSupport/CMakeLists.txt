configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.rc.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.rc")

include_external_msproject (GismoSupport
                           ${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport.vcxproj
                            )

set_target_properties(GismoSupport PROPERTIES
    FOLDER "Rhino THB plug-ins"  )							