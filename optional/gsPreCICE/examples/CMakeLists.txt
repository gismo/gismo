######################################################################
## CMakeLists.txt --- gsStructuralAnalysis/src
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Hugo Verhelst
## Copyright (C) 2023
######################################################################

# Collect all files in examples/
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# Collect all files in examples/subdir/
ALLSUBDIRS(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
FOREACH(subdir ${SUBDIRS})
  aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/${subdir} FILES)
  foreach(file ${FILES})
      add_gismo_executable(${file})
      get_filename_component(tarname ${file} NAME_WE) # name without extension
      set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
      set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
      # install the example executables (optionally)
      install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
  endforeach(file ${FILES})
ENDFOREACH()

# Collect all files in examples/subdir/subsubdir/
ALLSUBDIRS(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
FOREACH(subdir ${SUBDIRS})
  ALLSUBDIRS(SUBSUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
  FOREACH(subsubdir ${SUBSUBDIRS})
    aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/${subsubdir} FILES)
    foreach(file ${FILES})
        add_gismo_executable(${file})
        get_filename_component(tarname ${file} NAME_WE) # name without extension
        set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
        set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
        # install the example executables (optionally)
        install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
    endforeach(file ${FILES})
  ENDFOREACH()
ENDFOREACH()