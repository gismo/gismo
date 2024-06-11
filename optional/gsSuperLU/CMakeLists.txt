### CMakeLists.txt --- SuperLU library
## Author: Angelos Mantzaflaris
#
# The library is header-only, the extra code is downloaded in place of
# the source folder of G+Smo
#
######################################################################

## Spectra extension
project(gsSuperLUExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
#aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

include(gsFetch)
gismo_fetch_directory(SuperLUMT
  GIT_REPOSITORY https://github.com/xiaoyeli/superlu_mt.git
  #SVN_TRUST_CERT 1
  DESTINATION external
)

set(SuperLU_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/SuperLUMT-prefix/include/superlu_mt CACHE INTERNAL "")
set(SuperLU_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/SuperLUMT-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}superlu_mt_OPENMP${CMAKE_STATIC_LIBRARY_SUFFIX}
  CACHE INTERNAL "")

include(ExternalProject)
ExternalProject_Add(SuperLUMT
  BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/SuperLUMT
  SOURCE_DIR        ${gismo_externals}/SuperLUMT
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
             -DCMAKE_INSTALL_MESSAGE=NEVER
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -Denable_complex=OFF
             -Denable_complex16=OFF
             -Denable_examples=OFF
             -Denable_tests=OFF
             -Denable_single=OFF
             -DPLAT=_OPENMP
             # -DCMAKE_BUILD_TYPE=Release
  UPDATE_COMMAND     ""
  BUILD_BYPRODUCTS "${SuperLU_LIBRARY}"
  DOWNLOAD_COMMAND  ""
  #DOWNLOAD_NO_PROGRESS 0
  )

include_directories(SYSTEM ${SuperLU_INCLUDE_DIR})

set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${SuperLU_INCLUDE_DIR}
  CACHE INTERNAL "Gismo include directories")

set(gismo_LINKER ${gismo_LINKER} ${SuperLU_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")