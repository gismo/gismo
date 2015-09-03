 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2014 RICAM-Linz.
######################################################################

## Parasolid option
if(GISMO_WITH_PSOLID)

  add_gismo_executable(gsParasolid_test.cpp)
  mark_gismo_optional (gsParasolid_test.cpp)

  add_gismo_executable(gsParasolidBCURVE.cpp)
  mark_gismo_optional (gsParasolidBCURVE.cpp)

  add_gismo_executable(gsParasolidExportMesh.cpp)
  mark_gismo_optional (gsParasolidExportMesh.cpp)

  add_gismo_executable(gsParasolidExportTHB.cpp)
  mark_gismo_optional (gsParasolidExportTHB.cpp)
  
  #set_property(TARGET gsParasolid_test APPEND PROPERTY 
  #INCLUDE_DIRECTORIES ${PARASOLID_INCLUDE_DIR} )

endif(GISMO_WITH_PSOLID)


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
