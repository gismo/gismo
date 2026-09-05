######################################################################
## CMakeLists.txt --- gsCInterface
## This file is part of the G+Smo library.
##
## Type-blind core of the C interface: opaque handle types, the
## exception firewall (gsCError) and the generated umbrella
## <gsCInterface/Cgismo.h>. The C functions of a module live in
## src/<module>/cinterface/ and are compiled into the module when
## GISMO_BUILD_CINTERFACE is ON (see cmake/gsModule.cmake).
######################################################################

gismo_add_module(gsCInterface)

if(GISMO_BUILD_CINTERFACE)
  # plain-C examples double as the tests of the C interface
  aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples C_EXAMPLES)
  foreach(file ${C_EXAMPLES})
    if(file MATCHES "\\.c$")
      add_gismo_executable(${file})
      get_filename_component(tarname ${file} NAME_WE)
      set_property(TEST ${tarname} PROPERTY LABELS gsCInterface)
      set_target_properties(${tarname} PROPERTIES FOLDER "gsCInterface")
    endif()
  endforeach()
endif()
