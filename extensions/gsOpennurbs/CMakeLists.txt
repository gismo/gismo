### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2013 - RICAM-Linz.
######################################################################

## OpenNurbs extension
project(gsOpennurbsExtension)


# don't build shared libs 
set(BUILD_SHARED_LIBS OFF)

# build static libs by default
set(BUILD_STATIC_LIBS ON)

set(OPENNURBS_SOURCES
onurbs/opennurbs_3dm_attributes.cpp
onurbs/opennurbs_3dm_properties.cpp
onurbs/opennurbs_3dm_settings.cpp
onurbs/opennurbs_annotation2.cpp
onurbs/opennurbs_annotation.cpp
onurbs/opennurbs_arc.cpp
onurbs/opennurbs_arccurve.cpp
onurbs/opennurbs_archive.cpp
onurbs/opennurbs_array.cpp
onurbs/opennurbs_base32.cpp
onurbs/opennurbs_base64.cpp
onurbs/opennurbs_beam.cpp
onurbs/opennurbs_bezier.cpp
onurbs/opennurbs_beziervolume.cpp
onurbs/opennurbs_bitmap.cpp
onurbs/opennurbs_bounding_box.cpp
onurbs/opennurbs_box.cpp
onurbs/opennurbs_brep.cpp
onurbs/opennurbs_brep_extrude.cpp
onurbs/opennurbs_brep_io.cpp
onurbs/opennurbs_brep_isvalid.cpp
onurbs/opennurbs_brep_region.cpp
onurbs/opennurbs_brep_tools.cpp
onurbs/opennurbs_brep_v2valid.cpp
onurbs/opennurbs_circle.cpp
onurbs/opennurbs_color.cpp
onurbs/opennurbs_compress.cpp
onurbs/opennurbs_cone.cpp
onurbs/opennurbs_crc.cpp
onurbs/opennurbs_curve.cpp
onurbs/opennurbs_curveonsurface.cpp
onurbs/opennurbs_curveproxy.cpp
onurbs/opennurbs_cylinder.cpp
onurbs/opennurbs_defines.cpp
onurbs/opennurbs_detail.cpp
onurbs/opennurbs_dimstyle.cpp
onurbs/opennurbs_ellipse.cpp
onurbs/opennurbs_embedded_file.cpp
onurbs/opennurbs_error.cpp
onurbs/opennurbs_error_message.cpp
onurbs/opennurbs_evaluate_nurbs.cpp
onurbs/opennurbs_extensions.cpp
onurbs/opennurbs_font.cpp
onurbs/opennurbs_fsp.cpp
onurbs/opennurbs_geometry.cpp
onurbs/opennurbs_group.cpp
onurbs/opennurbs_hatch.cpp
onurbs/opennurbs_instance.cpp
onurbs/opennurbs_intersect.cpp
onurbs/opennurbs_knot.cpp
onurbs/opennurbs_layer.cpp
onurbs/opennurbs_light.cpp
onurbs/opennurbs_line.cpp
onurbs/opennurbs_linecurve.cpp
onurbs/opennurbs_linetype.cpp
onurbs/opennurbs_lookup.cpp
onurbs/opennurbs_material.cpp
onurbs/opennurbs_math.cpp
onurbs/opennurbs_matrix.cpp
onurbs/opennurbs_mesh.cpp
onurbs/opennurbs_mesh_ngon.cpp
onurbs/opennurbs_mesh_tools.cpp
onurbs/opennurbs_morph.cpp
onurbs/opennurbs_nurbscurve.cpp
onurbs/opennurbs_nurbssurface.cpp
onurbs/opennurbs_nurbsvolume.cpp
onurbs/opennurbs_object.cpp
onurbs/opennurbs_object_history.cpp
onurbs/opennurbs_objref.cpp
onurbs/opennurbs_offsetsurface.cpp
onurbs/opennurbs_optimize.cpp
onurbs/opennurbs_plane.cpp
onurbs/opennurbs_planesurface.cpp
onurbs/opennurbs_pluginlist.cpp
onurbs/opennurbs_pointcloud.cpp
onurbs/opennurbs_point.cpp
onurbs/opennurbs_pointgeometry.cpp
onurbs/opennurbs_pointgrid.cpp
onurbs/opennurbs_polycurve.cpp
onurbs/opennurbs_polyedgecurve.cpp
onurbs/opennurbs_polyline.cpp
onurbs/opennurbs_polylinecurve.cpp
onurbs/opennurbs_rand.cpp
onurbs/opennurbs_revsurface.cpp
onurbs/opennurbs_rtree.cpp
onurbs/opennurbs_sort.cpp
onurbs/opennurbs_sphere.cpp
onurbs/opennurbs_string.cpp
onurbs/opennurbs_sum.cpp
onurbs/opennurbs_sumsurface.cpp
onurbs/opennurbs_surface.cpp
onurbs/opennurbs_surfaceproxy.cpp
onurbs/opennurbs_textlog.cpp
onurbs/opennurbs_torus.cpp
onurbs/opennurbs_unicode.cpp
onurbs/opennurbs_userdata.cpp
onurbs/opennurbs_uuid.cpp
onurbs/opennurbs_viewport.cpp
onurbs/opennurbs_workspace.cpp
onurbs/opennurbs_wstring.cpp
onurbs/opennurbs_xform.cpp
onurbs/opennurbs_zlib.cpp
onurbs/opennurbs_zlib_memory.cpp
onurbs/opennurbs_memory.c
onurbs/opennurbs_memory_util.c
# zlib
onurbs/zlib/adler32.c
onurbs/zlib/compress.c
onurbs/zlib/crc32.c
onurbs/zlib/deflate.c
onurbs/zlib/infback.c
onurbs/zlib/inffast.c
onurbs/zlib/inflate.c
onurbs/zlib/inftrees.c
onurbs/zlib/trees.c
onurbs/zlib/uncompr.c
onurbs/zlib/zutil.c
# 
gsReadOpenNurbs.cpp
  )

if (BUILD_SHARED_LIBS)
set(OPENNURBS_SOURCES ${OPENNURBS_SOURCES} onurbs/opennurbs_dll.cpp)
endif()

set(OPENNURBS_SOURCES_UNUSED
onurbs/opennurbs_gl.cpp
  )


SET(OPENNURBS_HEADERS
examples_linking_pragmas.h
onurbs/opennurbs_3dm_attributes.h
onurbs/opennurbs_3dm.h
onurbs/opennurbs_3dm_properties.h
onurbs/opennurbs_3dm_settings.h
onurbs/opennurbs_annotation2.h
onurbs/opennurbs_annotation.h
onurbs/opennurbs_arccurve.h
onurbs/opennurbs_arc.h
onurbs/opennurbs_archive.h
onurbs/opennurbs_array_defs.h
onurbs/opennurbs_array.h
onurbs/opennurbs_base32.h
onurbs/opennurbs_base64.h
onurbs/opennurbs_beam.h
onurbs/opennurbs_bezier.h
onurbs/opennurbs_bitmap.h
onurbs/opennurbs_bounding_box.h
onurbs/opennurbs_box.h
onurbs/opennurbs_brep.h
onurbs/opennurbs_circle.h
onurbs/opennurbs_color.h
onurbs/opennurbs_compress.h
onurbs/opennurbs_cone.h
onurbs/opennurbs_crc.h
onurbs/opennurbs_curve.h
onurbs/opennurbs_curveonsurface.h
onurbs/opennurbs_curveproxy.h
onurbs/opennurbs_cylinder.h
onurbs/opennurbs_defines.h
onurbs/opennurbs_detail.h
onurbs/opennurbs_dimstyle.h
onurbs/opennurbs_ellipse.h
onurbs/opennurbs_error.h
onurbs/opennurbs_evaluate_nurbs.h
onurbs/opennurbs_extensions.h
onurbs/opennurbs_font.h
onurbs/opennurbs_fpoint.h
onurbs/opennurbs_fsp_defs.h
onurbs/opennurbs_fsp.h
onurbs/opennurbs_geometry.h
onurbs/opennurbs_group.h
opennurbs.h
onurbs/opennurbs_hatch.h
onurbs/opennurbs_hsort_template.h
onurbs/opennurbs_instance.h
onurbs/opennurbs_intersect.h
onurbs/opennurbs_knot.h
onurbs/opennurbs_layer.h
onurbs/opennurbs_light.h
onurbs/opennurbs_linecurve.h
onurbs/opennurbs_line.h
onurbs/opennurbs_linestyle.h
onurbs/opennurbs_linetype.h
onurbs/opennurbs_lookup.h
onurbs/opennurbs_mapchan.h
onurbs/opennurbs_material.h
onurbs/opennurbs_math.h
onurbs/opennurbs_matrix.h
onurbs/opennurbs_memory.h
onurbs/opennurbs_mesh.h
onurbs/opennurbs_nurbscurve.h
onurbs/opennurbs_nurbssurface.h
onurbs/opennurbs_object.h
onurbs/opennurbs_object_history.h
onurbs/opennurbs_objref.h
onurbs/opennurbs_offsetsurface.h
onurbs/opennurbs_optimize.h
onurbs/opennurbs_plane.h
onurbs/opennurbs_planesurface.h
onurbs/opennurbs_pluginlist.h
onurbs/opennurbs_pointcloud.h
onurbs/opennurbs_pointgeometry.h
onurbs/opennurbs_pointgrid.h
onurbs/opennurbs_point.h
onurbs/opennurbs_polycurve.h
onurbs/opennurbs_polyedgecurve.h
onurbs/opennurbs_polylinecurve.h
onurbs/opennurbs_polyline.h
onurbs/opennurbs_qsort_template.h
onurbs/opennurbs_rand.h
onurbs/opennurbs_rendering.h
onurbs/opennurbs_revsurface.h
onurbs/opennurbs_rtree.h
onurbs/opennurbs_sphere.h
onurbs/opennurbs_string.h
onurbs/opennurbs_sumsurface.h
onurbs/opennurbs_surface.h
onurbs/opennurbs_surfaceproxy.h
onurbs/opennurbs_system.h
onurbs/opennurbs_textlog.h
onurbs/opennurbs_texture.h
onurbs/opennurbs_texture_mapping.h
onurbs/opennurbs_torus.h
onurbs/opennurbs_unicode.h
onurbs/opennurbs_userdata.h
onurbs/opennurbs_uuid.h
onurbs/opennurbs_version.h
onurbs/opennurbs_viewport.h
onurbs/opennurbs_workspace.h
onurbs/opennurbs_xform.h
onurbs/opennurbs_zlib.h
# zlib
onurbs/zlib/crc32.h
onurbs/zlib/deflate.h
onurbs/zlib/inffast.h
onurbs/zlib/inffixed.h
onurbs/zlib/inflate.h
onurbs/zlib/inftrees.h
onurbs/zlib/trees.h
onurbs/zlib/zconf.h
onurbs/zlib/zlib.h
onurbs/zlib/zutil.h
# 
gsReadOpenNurbs.h
  )

set(OPENNURBS_UNUSED
onurbs/opennurbs_gl.h
onurbs/opennurbs_precompiledheader.cpp
  )


SET(ON_DLL_HEADERS
  onurbs/opennurbs_dll_resource.h
  )

INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES(emmintrin.h HAVE_EMMINTRIN_H)
CHECK_INCLUDE_FILES(ieeefp.h HAVE_IEEEFP_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
if(HAVE_EMMINTRIN_H)
  add_definitions(-DHAVE_EMMINTRIN_H)
endif(HAVE_EMMINTRIN_H)
if(HAVE_IEEEFP_H)
  add_definitions(-DHAVE_IEEEFP_H)
endif(HAVE_IEEEFP_H)
if(HAVE_SYS_STAT_H)
  add_definitions(-DHAVE_SYS_STAT_H)
endif(HAVE_SYS_STAT_H)

# Ideally, we should only need ZLIB here
#include_directories(
#  ${${CMAKE_PROJECT_NAME}_SOURCE_DIR}/onurbs/zlib
#  )

# Apply same configuration as G+Smo
include( gsConfig)

#if (MSVC)
#  add_definitions(
#    -DON_COMPILING_OPENNURBS
#    -DON_DLL_EXPORTS
#    -DON_DLL_IMPORTS
#    )
#  SET(OPENNURBS_HEADERS "${OPENNURBS_HEADERS};${ON_DLL_HEADERS}")
#else (MSVC)
add_definitions(-DON_COMPILING_OPENNURBS)
#endif(MSVC)


if(MSVC)
  add_definitions(-DUNICODE)
  set(OPENNURBS_LINKLIBRARIES ${OPENNURBS_LINKLIBRARIES} Rpcrt4 ) # needed ?
endif(MSVC)

# Do not re-compile G+Smo in this unit
if ( NOT GISMO_BUILD_LIB)
   add_definitions(-DGISMO_BUILD_LIB)
endif()

IF(BUILD_STATIC_LIBS)
    add_library(${PROJECT_NAME} OBJECT ${OPENNURBS_SOURCES}) #STATIC
    #target_link_libraries(${PROJECT_NAME} ${OPENNURBS_LINKLIBRARIES}) # needed ?
    IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
        SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS -fPIC)
    ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
    IF(NOT WIN32)
      SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}")
    ENDIF(NOT WIN32)
ELSEIF(BUILD_SHARED_LIBS)
  add_library(${PROJECT_NAME} SHARED ${OPENNURBS_SOURCES})
  #target_link_libraries(${PROJECT_NAME} ${OPENNURBS_LINKLIBRARIES}) # needed ?
  set_property(TARGET ${PROJECT_NAME} PROPERTY LINK_INTERFACE_LIBRARIES "")
  SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION "2010.04.95")
ENDIF(BUILD_STATIC_LIBS)

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:gsOpennurbsExtension>
    CACHE INTERNAL "Gismo extensions to be included")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsOpennurbs/
        FILES_MATCHING PATTERN "*.h")
