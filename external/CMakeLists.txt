######################################################################
### CMakeLists.txt --- external
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

set(${PROJECT_NAME}_externals "${CMAKE_CURRENT_SOURCE_DIR}"
    CACHE INTERNAL "Externals directory")

include(gszlib.cmake)
include(gszstream.cmake)

if(GISMO_BUILD_CPPLOT)
#add_subdirectory(boost)
add_subdirectory(cpplot)
endif(GISMO_BUILD_CPPLOT)
