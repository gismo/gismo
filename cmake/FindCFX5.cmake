
if (CFX5_INCLUDES AND CFX5_LIBRARIES)
  set(CFX5_FIND_QUIETLY TRUE)
endif (CFX5_INCLUDES AND CFX5_LIBRARIES)

find_path(CFX5_INCLUDES
  NAMES
  cfxids.h
  PATHS
  $ENV{CFX5ROOT}/include
  ${INCLUDE_INSTALL_DIR}
)


if (WIN32)
  find_library(CFX5_MESHEXPORT_LIBRARY meshexport PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_RATLAS_API_LIBRARY ratlas_api PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_RATLAS_LIBRARY     ratlas     PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_PGTAPI_LIBRARY     pgtapi     PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_UNITS_LIBRARY      units      PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_CCLAPILT_LIBRARY   cclapilt   PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_IO_LIBRARY         io         PATHS $ENV{CFX5ROOT}/lib/winnt-amd64 ${LIB_INSTALL_DIR})
elseif(UNIX AND NOT APPLE)
  find_library(CFX5_MESHEXPORT_LIBRARY meshexport PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_RATLAS_API_LIBRARY ratlas_api PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_RATLAS_LIBRARY     ratlas     PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_PGTAPI_LIBRARY     pgtapi     PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_UNITS_LIBRARY      units      PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_CCLAPILT_LIBRARY   cclapilt   PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
  find_library(CFX5_IO_LIBRARY         io         PATHS $ENV{CFX5ROOT}/lib/linux-amd64 ${LIB_INSTALL_DIR})
else()
  error("CFX5-support is not available for your platform.")
endif()

set(CFX5_LIBRARIES ${CFX5_MESHEXPORT_LIBRARY} ${CFX5_RATLAS_API_LIBRARY} ${CFX5_RATLAS_LIBRARY} ${CFX5_PGTAPI_LIBRARY} ${CFX5_UNITS_LIBRARY} ${CFX5_CCLAPILT_LIBRARY} ${CFX5_IO_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFX5 DEFAULT_MSG
                                  CFX5_INCLUDES CFX5_LIBRARIES)

mark_as_advanced(CFX5_INCLUDES CFX5_LIBRARIES)
