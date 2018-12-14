# Fetching and downloading gmp

gismo_fetch_directory(gmp
  URL https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2
  DESTINATION external)

include(ExternalProject)

ExternalProject_Add(gmp_build
  SOURCE_DIR ${gismo_SOURCE_DIR}/external/gmp
  BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}/gmp_build
  CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${gismo_SOURCE_DIR}/external/gmp/configure --enable-cxx --prefix=<INSTALL_DIR> )

ExternalProject_Get_Property(gmp_build install_dir)
set(GS_GMP_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
set(GS_GMP_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
set(GS_GMP_LIBRARY ${GS_GMP_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}gmp${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE INTERNAL "")

# Ensure this compiles prior to gismo default modules
add_dependencies(gsDepends gmp_build)

set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GS_GMP_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
set(gismo_LINKER ${gismo_LINKER} ${GS_GMP_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

# to do
#install( FILES ${HEADERS} DESTINATION include/ OPTIONAL)
