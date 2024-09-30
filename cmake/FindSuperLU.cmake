
if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
  set(SUPERLU_FIND_QUIETLY TRUE)
endif (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

find_package(BLAS QUIET)

if(BLAS_FOUND)

  # Search in user-defined paths only
  find_path(SUPERLU_INCLUDES
    NAMES "supermatrix.h"
    PATHS
    $ENV{SUPERLUDIR}
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES "superlu" "SuperLU" "include/superlu" "include" "SRC"
    NO_DEFAULT_PATH
    )

  # Search in default paths
  find_path(SUPERLU_INCLUDES
    NAMES "supermatrix.h"
    PATH_SUFFIXES "superlu" "SuperLU" "include/superlu" "include" "SRC"
    )

  # Search in user-defined paths only
  find_library(SUPERLU_LIBRARIES
    NAMES "superlu_4.3" "superlu_4.2" "superlu_4.1" "superlu_4.0" "superlu"
    PATHS
    $ENV{SUPERLUDIR}
    ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib32" "lib64"
    NO_DEFAULT_PATH
    )

  # Search in default paths
  find_library(SUPERLU_LIBRARIES
    NAMES "superlu_4.3" "superlu_4.2" "superlu_4.1" "superlu_4.0" "superlu"
    PATHS
    $ENV{SUPERLUDIR}
    ${LIB_INSTALL_DIR}
    PATH_SUFFIXES "lib" "lib32" "lib64"
    )
  
  if(SUPERLU_LIBRARIES)
    set(SUPERLU_LIBRARIES ${SUPERLU_LIBRARIES} ${BLAS_LIBRARIES})
  endif(SUPERLU_LIBRARIES)
  
endif(BLAS_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUPERLU DEFAULT_MSG
  SUPERLU_INCLUDES SUPERLU_LIBRARIES)

mark_as_advanced(SUPERLU_INCLUDES SUPERLU_LIBRARIES)
