######################################################################
### CMakeLists.txt --- gsOpennurbs
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

## OpenNurbs extension
project(gsOpennurbs)

set(OPENNURBS_SOURCES
${gismo_externals}/onurbs/opennurbs_3dm_attributes.cpp
${gismo_externals}/onurbs/opennurbs_3dm_properties.cpp
${gismo_externals}/onurbs/opennurbs_3dm_settings.cpp
${gismo_externals}/onurbs/opennurbs_annotation2.cpp
${gismo_externals}/onurbs/opennurbs_annotation.cpp
${gismo_externals}/onurbs/opennurbs_arc.cpp
${gismo_externals}/onurbs/opennurbs_arccurve.cpp
${gismo_externals}/onurbs/opennurbs_archive.cpp
${gismo_externals}/onurbs/opennurbs_array.cpp
${gismo_externals}/onurbs/opennurbs_base32.cpp
${gismo_externals}/onurbs/opennurbs_base64.cpp
${gismo_externals}/onurbs/opennurbs_beam.cpp
${gismo_externals}/onurbs/opennurbs_bezier.cpp
${gismo_externals}/onurbs/opennurbs_beziervolume.cpp
${gismo_externals}/onurbs/opennurbs_bitmap.cpp
${gismo_externals}/onurbs/opennurbs_bounding_box.cpp
${gismo_externals}/onurbs/opennurbs_box.cpp
${gismo_externals}/onurbs/opennurbs_brep.cpp
${gismo_externals}/onurbs/opennurbs_brep_extrude.cpp
${gismo_externals}/onurbs/opennurbs_brep_io.cpp
${gismo_externals}/onurbs/opennurbs_brep_isvalid.cpp
${gismo_externals}/onurbs/opennurbs_brep_region.cpp
${gismo_externals}/onurbs/opennurbs_brep_tools.cpp
${gismo_externals}/onurbs/opennurbs_brep_v2valid.cpp
${gismo_externals}/onurbs/opennurbs_circle.cpp
${gismo_externals}/onurbs/opennurbs_color.cpp
${gismo_externals}/onurbs/opennurbs_compress.cpp
${gismo_externals}/onurbs/opennurbs_cone.cpp
${gismo_externals}/onurbs/opennurbs_crc.cpp
${gismo_externals}/onurbs/opennurbs_curve.cpp
${gismo_externals}/onurbs/opennurbs_curveonsurface.cpp
${gismo_externals}/onurbs/opennurbs_curveproxy.cpp
${gismo_externals}/onurbs/opennurbs_cylinder.cpp
${gismo_externals}/onurbs/opennurbs_defines.cpp
${gismo_externals}/onurbs/opennurbs_detail.cpp
${gismo_externals}/onurbs/opennurbs_dimstyle.cpp
${gismo_externals}/onurbs/opennurbs_ellipse.cpp
${gismo_externals}/onurbs/opennurbs_embedded_file.cpp
${gismo_externals}/onurbs/opennurbs_error.cpp
${gismo_externals}/onurbs/opennurbs_error_message.cpp
${gismo_externals}/onurbs/opennurbs_evaluate_nurbs.cpp
${gismo_externals}/onurbs/opennurbs_extensions.cpp
${gismo_externals}/onurbs/opennurbs_font.cpp
${gismo_externals}/onurbs/opennurbs_fsp.cpp
${gismo_externals}/onurbs/opennurbs_geometry.cpp
${gismo_externals}/onurbs/opennurbs_group.cpp
${gismo_externals}/onurbs/opennurbs_hatch.cpp
${gismo_externals}/onurbs/opennurbs_instance.cpp
${gismo_externals}/onurbs/opennurbs_intersect.cpp
${gismo_externals}/onurbs/opennurbs_knot.cpp
${gismo_externals}/onurbs/opennurbs_layer.cpp
${gismo_externals}/onurbs/opennurbs_light.cpp
${gismo_externals}/onurbs/opennurbs_line.cpp
${gismo_externals}/onurbs/opennurbs_linecurve.cpp
${gismo_externals}/onurbs/opennurbs_linetype.cpp
${gismo_externals}/onurbs/opennurbs_lookup.cpp
${gismo_externals}/onurbs/opennurbs_material.cpp
${gismo_externals}/onurbs/opennurbs_math.cpp
${gismo_externals}/onurbs/opennurbs_matrix.cpp
${gismo_externals}/onurbs/opennurbs_mesh.cpp
${gismo_externals}/onurbs/opennurbs_mesh_ngon.cpp
${gismo_externals}/onurbs/opennurbs_mesh_tools.cpp
${gismo_externals}/onurbs/opennurbs_morph.cpp
${gismo_externals}/onurbs/opennurbs_nurbscurve.cpp
${gismo_externals}/onurbs/opennurbs_nurbssurface.cpp
${gismo_externals}/onurbs/opennurbs_nurbsvolume.cpp
${gismo_externals}/onurbs/opennurbs_object.cpp
${gismo_externals}/onurbs/opennurbs_object_history.cpp
${gismo_externals}/onurbs/opennurbs_objref.cpp
${gismo_externals}/onurbs/opennurbs_offsetsurface.cpp
${gismo_externals}/onurbs/opennurbs_optimize.cpp
${gismo_externals}/onurbs/opennurbs_plane.cpp
${gismo_externals}/onurbs/opennurbs_planesurface.cpp
${gismo_externals}/onurbs/opennurbs_pluginlist.cpp
${gismo_externals}/onurbs/opennurbs_pointcloud.cpp
${gismo_externals}/onurbs/opennurbs_point.cpp
${gismo_externals}/onurbs/opennurbs_pointgeometry.cpp
${gismo_externals}/onurbs/opennurbs_pointgrid.cpp
${gismo_externals}/onurbs/opennurbs_polycurve.cpp
${gismo_externals}/onurbs/opennurbs_polyedgecurve.cpp
${gismo_externals}/onurbs/opennurbs_polyline.cpp
${gismo_externals}/onurbs/opennurbs_polylinecurve.cpp
${gismo_externals}/onurbs/opennurbs_rand.cpp
${gismo_externals}/onurbs/opennurbs_revsurface.cpp
${gismo_externals}/onurbs/opennurbs_rtree.cpp
${gismo_externals}/onurbs/opennurbs_sort.cpp
${gismo_externals}/onurbs/opennurbs_sphere.cpp
${gismo_externals}/onurbs/opennurbs_string.cpp
${gismo_externals}/onurbs/opennurbs_sum.cpp
${gismo_externals}/onurbs/opennurbs_sumsurface.cpp
${gismo_externals}/onurbs/opennurbs_surface.cpp
${gismo_externals}/onurbs/opennurbs_surfaceproxy.cpp
${gismo_externals}/onurbs/opennurbs_textlog.cpp
${gismo_externals}/onurbs/opennurbs_torus.cpp
${gismo_externals}/onurbs/opennurbs_unicode.cpp
${gismo_externals}/onurbs/opennurbs_userdata.cpp
${gismo_externals}/onurbs/opennurbs_uuid.cpp
${gismo_externals}/onurbs/opennurbs_viewport.cpp
${gismo_externals}/onurbs/opennurbs_workspace.cpp
${gismo_externals}/onurbs/opennurbs_wstring.cpp
${gismo_externals}/onurbs/opennurbs_xform.cpp
${gismo_externals}/onurbs/opennurbs_zlib.cpp
${gismo_externals}/onurbs/opennurbs_zlib_memory.cpp
${gismo_externals}/onurbs/opennurbs_quaternion.cpp
${gismo_externals}/onurbs/opennurbs_photogrammetry.cpp
${gismo_externals}/onurbs/opennurbs_memory.c
${gismo_externals}/onurbs/opennurbs_memory_util.c
# zlib
#${gismo_externals}/onurbs/zlib/adler32.c
#${gismo_externals}/onurbs/zlib/compress.c
#${gismo_externals}/onurbs/zlib/crc32.c
#${gismo_externals}/onurbs/zlib/deflate.c
#${gismo_externals}/onurbs/zlib/infback.c
#${gismo_externals}/onurbs/zlib/inffast.c
#${gismo_externals}/onurbs/zlib/inflate.c
#${gismo_externals}/onurbs/zlib/inftrees.c
#${gismo_externals}/onurbs/zlib/trees.c
#${gismo_externals}/onurbs/zlib/uncompr.c
#${gismo_externals}/onurbs/zlib/zutil.c
# 
gsReadOpenNurbs.cpp
#gsWriteOpenNurbs_.cpp
)

if (GISMO_BUILD_LIB)
set(OPENNURBS_SOURCES ${OPENNURBS_SOURCES} gsWriteOpenNurbs_.cpp)
endif()

#if (BUILD_SHARED_LIBS)
#set(OPENNURBS_SOURCES ${OPENNURBS_SOURCES} ${gismo_externals}/onurbs/opennurbs_dll.cpp)
#endif()
#set(OPENNURBS_SOURCES_UNUSED
#${gismo_externals}/onurbs/opennurbs_gl.cpp
#  )

SET(OPENNURBS_HEADERS
examples_linking_pragmas.h
${gismo_externals}/onurbs/opennurbs_3dm_attributes.h
${gismo_externals}/onurbs/opennurbs_3dm.h
${gismo_externals}/onurbs/opennurbs_3dm_properties.h
${gismo_externals}/onurbs/opennurbs_3dm_settings.h
${gismo_externals}/onurbs/opennurbs_annotation2.h
${gismo_externals}/onurbs/opennurbs_annotation.h
${gismo_externals}/onurbs/opennurbs_arccurve.h
${gismo_externals}/onurbs/opennurbs_arc.h
${gismo_externals}/onurbs/opennurbs_archive.h
${gismo_externals}/onurbs/opennurbs_array_defs.h
${gismo_externals}/onurbs/opennurbs_array.h
${gismo_externals}/onurbs/opennurbs_base32.h
${gismo_externals}/onurbs/opennurbs_base64.h
${gismo_externals}/onurbs/opennurbs_beam.h
${gismo_externals}/onurbs/opennurbs_bezier.h
${gismo_externals}/onurbs/opennurbs_bitmap.h
${gismo_externals}/onurbs/opennurbs_bounding_box.h
${gismo_externals}/onurbs/opennurbs_box.h
${gismo_externals}/onurbs/opennurbs_brep.h
${gismo_externals}/onurbs/opennurbs_circle.h
${gismo_externals}/onurbs/opennurbs_color.h
${gismo_externals}/onurbs/opennurbs_compress.h
${gismo_externals}/onurbs/opennurbs_cone.h
${gismo_externals}/onurbs/opennurbs_crc.h
${gismo_externals}/onurbs/opennurbs_curve.h
${gismo_externals}/onurbs/opennurbs_curveonsurface.h
${gismo_externals}/onurbs/opennurbs_curveproxy.h
${gismo_externals}/onurbs/opennurbs_cylinder.h
${gismo_externals}/onurbs/opennurbs_defines.h
${gismo_externals}/onurbs/opennurbs_detail.h
${gismo_externals}/onurbs/opennurbs_dimstyle.h
${gismo_externals}/onurbs/opennurbs_ellipse.h
${gismo_externals}/onurbs/opennurbs_error.h
${gismo_externals}/onurbs/opennurbs_evaluate_nurbs.h
${gismo_externals}/onurbs/opennurbs_extensions.h
${gismo_externals}/onurbs/opennurbs_font.h
${gismo_externals}/onurbs/opennurbs_fpoint.h
${gismo_externals}/onurbs/opennurbs_fsp_defs.h
${gismo_externals}/onurbs/opennurbs_fsp.h
${gismo_externals}/onurbs/opennurbs_geometry.h
${gismo_externals}/onurbs/opennurbs_group.h
opennurbs.h
${gismo_externals}/onurbs/opennurbs_hatch.h
${gismo_externals}/onurbs/opennurbs_hsort_template.h
${gismo_externals}/onurbs/opennurbs_instance.h
${gismo_externals}/onurbs/opennurbs_intersect.h
${gismo_externals}/onurbs/opennurbs_knot.h
${gismo_externals}/onurbs/opennurbs_layer.h
${gismo_externals}/onurbs/opennurbs_light.h
${gismo_externals}/onurbs/opennurbs_linecurve.h
${gismo_externals}/onurbs/opennurbs_line.h
${gismo_externals}/onurbs/opennurbs_linestyle.h
${gismo_externals}/onurbs/opennurbs_linetype.h
${gismo_externals}/onurbs/opennurbs_lookup.h
${gismo_externals}/onurbs/opennurbs_mapchan.h
${gismo_externals}/onurbs/opennurbs_material.h
${gismo_externals}/onurbs/opennurbs_math.h
${gismo_externals}/onurbs/opennurbs_matrix.h
${gismo_externals}/onurbs/opennurbs_memory.h
${gismo_externals}/onurbs/opennurbs_mesh.h
${gismo_externals}/onurbs/opennurbs_nurbscurve.h
${gismo_externals}/onurbs/opennurbs_nurbssurface.h
${gismo_externals}/onurbs/opennurbs_object.h
${gismo_externals}/onurbs/opennurbs_object_history.h
${gismo_externals}/onurbs/opennurbs_objref.h
${gismo_externals}/onurbs/opennurbs_offsetsurface.h
${gismo_externals}/onurbs/opennurbs_optimize.h
${gismo_externals}/onurbs/opennurbs_plane.h
${gismo_externals}/onurbs/opennurbs_planesurface.h
${gismo_externals}/onurbs/opennurbs_pluginlist.h
${gismo_externals}/onurbs/opennurbs_pointcloud.h
${gismo_externals}/onurbs/opennurbs_pointgeometry.h
${gismo_externals}/onurbs/opennurbs_pointgrid.h
${gismo_externals}/onurbs/opennurbs_point.h
${gismo_externals}/onurbs/opennurbs_polycurve.h
${gismo_externals}/onurbs/opennurbs_polyedgecurve.h
${gismo_externals}/onurbs/opennurbs_polylinecurve.h
${gismo_externals}/onurbs/opennurbs_polyline.h
${gismo_externals}/onurbs/opennurbs_qsort_template.h
${gismo_externals}/onurbs/opennurbs_rand.h
${gismo_externals}/onurbs/opennurbs_rendering.h
${gismo_externals}/onurbs/opennurbs_revsurface.h
${gismo_externals}/onurbs/opennurbs_rtree.h
${gismo_externals}/onurbs/opennurbs_sphere.h
${gismo_externals}/onurbs/opennurbs_string.h
${gismo_externals}/onurbs/opennurbs_sumsurface.h
${gismo_externals}/onurbs/opennurbs_surface.h
${gismo_externals}/onurbs/opennurbs_surfaceproxy.h
${gismo_externals}/onurbs/opennurbs_system.h
${gismo_externals}/onurbs/opennurbs_textlog.h
${gismo_externals}/onurbs/opennurbs_texture.h
${gismo_externals}/onurbs/opennurbs_texture_mapping.h
${gismo_externals}/onurbs/opennurbs_torus.h
${gismo_externals}/onurbs/opennurbs_unicode.h
${gismo_externals}/onurbs/opennurbs_userdata.h
${gismo_externals}/onurbs/opennurbs_uuid.h
${gismo_externals}/onurbs/opennurbs_version.h
${gismo_externals}/onurbs/opennurbs_viewport.h
${gismo_externals}/onurbs/opennurbs_workspace.h
${gismo_externals}/onurbs/opennurbs_xform.h
${gismo_externals}/onurbs/opennurbs_zlib.h
#${gismo_externals}/onurbs/opennurbs_quaternion.h
${gismo_externals}/onurbs/opennurbs_photogrammetry.h
# zlib
#${gismo_externals}/onurbs/zlib/crc32.h
#${gismo_externals}/onurbs/zlib/deflate.h
#${gismo_externals}/onurbs/zlib/inffast.h
#${gismo_externals}/onurbs/zlib/inffixed.h
#${gismo_externals}/onurbs/zlib/inflate.h
#${gismo_externals}/onurbs/zlib/inftrees.h
#${gismo_externals}/onurbs/zlib/trees.h
#${gismo_externals}/onurbs/zlib/zconf.h
#${gismo_externals}/onurbs/zlib/zlib.h
#${gismo_externals}/onurbs/zlib/zutil.h
# 
gsReadOpenNurbs.h
gsWriteOpenNurbs.h
)

set(OPENNURBS_UNUSED
${gismo_externals}/onurbs/opennurbs_gl.h
${gismo_externals}/onurbs/opennurbs_precompiledheader.cpp
  )


SET(ON_DLL_HEADERS
  ${gismo_externals}/onurbs/opennurbs_dll_resource.h
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

# Apply same configuration as G+Smo
include(gsConfig)

#if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
#  add_definitions(
#    -DON_COMPILING_OPENNURBS
#    -DON_DLL_EXPORTS
#    -DON_DLL_IMPORTS
#    )
#  SET(OPENNURBS_HEADERS "${OPENNURBS_HEADERS};${ON_DLL_HEADERS}")
#else ()
add_definitions(-DON_COMPILING_OPENNURBS)
#endif()


if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  add_definitions(-DUNICODE)
  add_definitions(
    -wd4458 # Suppress 'declaration of '..' hides class member'	~2400 warns
    -wd4456 # declaration of '..' hides previous local declaration ~37 warns
    )
else()
  add_definitions(-DON_NO_WINDOWS) #for Mingw
endif()

if (GISMO_ZLIB_STATIC)
  #message("Using prefixed Zlib")
  add_definitions(-DZ_PREFIX)   # use prefixed zlib
  add_definitions(-DMY_ZCALLOC) # zlib uses oncalloc() and onfree()
endif()

add_library(${PROJECT_NAME} OBJECT ${OPENNURBS_SOURCES}) #STATIC
#target_link_libraries(${PROJECT_NAME} Rpcrt4) # Note: linked with gismolib instead

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/W0")
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-w")
endif()
#In CMake 2.8.12 the above can be:
#target_compile_options(${PROJECT_NAME} PRIVATE $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-w> $<$<CXX_COMPILER_ID:MSVC>: /W0>)

set_target_properties(${PROJECT_NAME} 
PROPERTIES POSITION_INDEPENDENT_CODE ON
COMPILE_DEFINITIONS gismo_EXPORTS
FOLDER "G+Smo extensions")

#SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION "2010.04.95")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")
if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" )
set(gismo_LINKER ${gismo_LINKER} Rpcrt4
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "*.h")
