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
  
  #set_property(TARGET gsParasolid_test APPEND PROPERTY 
  #INCLUDE_DIRECTORIES ${PARASOLID_INCLUDE_DIR} )

endif(GISMO_WITH_PSOLID)

if(GISMO_WITH_SUPERLU)
  add_gismo_executable(gsSuperLU_test.cpp)
  mark_gismo_optional (gsSuperLU_test.cpp)
endif(GISMO_WITH_SUPERLU)

if(GISMO_WITH_PARDISO)
  add_gismo_executable(gsPardiso_test.cpp)
  mark_gismo_optional (gsPardiso_test.cpp)
endif(GISMO_WITH_PARDISO)


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
