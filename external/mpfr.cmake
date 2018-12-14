# Fetching and downloading mpfr

gismo_fetch_directory(mpfr
  URL https://www.mpfr.org/mpfr-current/mpfr-4.0.1.tar.gz
  DESTINATION external)

include(ExternalProject)

ExternalProject_Add(mpfr_build
  SOURCE_DIR ${gismo_SOURCE_DIR}/external/mpfr
  BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}/mpfr_build
  CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${gismo_SOURCE_DIR}/external/mpfr/configure --enable-cxx --prefix=<INSTALL_DIR> )

ExternalProject_Get_Property(mpfr_build install_dir)
set(GS_MPFR_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
set(GS_MPFR_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
set(GS_MPFR_LIBRARY ${GS_MPFR_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}mpfr${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE INTERNAL "")

# Ensure this compiles prior to gismo default modules
add_dependencies(gsDepends mpfr_build)

set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GS_MPFR_INCLUDE_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include directories")
set(gismo_LINKER ${gismo_LINKER} ${GS_MPFR_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

# to do
#install( FILES ${HEADERS} DESTINATION include/ OPTIONAL)
