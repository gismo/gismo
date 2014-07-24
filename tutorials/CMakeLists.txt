 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

# Collect source file names
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)

if(NOT GISMO_BUILD_CPPLOT)
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/gsPlot_test.cpp)
endif(NOT GISMO_BUILD_CPPLOT)

foreach(file ${FILES})
    add_gismo_test(${file})
endforeach(file ${FILES})

if(GISMO_BUILD_CPPLOT)
# Link with cpplot for plotting
target_link_libraries(gsPlot_test 
cpplot
#boostpart
)
endif(GISMO_BUILD_CPPLOT)

if(GISMO_BUILD_QT_APP)
   	message(STATUS "found QtGismo")
	add_subdirectory( QtGismo)
endif(GISMO_BUILD_QT_APP)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

