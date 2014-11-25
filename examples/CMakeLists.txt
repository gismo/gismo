### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

# Collect sourcer file names
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
#message("found files: ${FILES}")

foreach(file ${FILES})
    add_gismo_example(${file})
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

