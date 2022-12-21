######################################################################
## gsOptions.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

message ("Configuration (cmake ${CMAKE_VERSION}):")

message ("  Source folder:          ${CMAKE_SOURCE_DIR}")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  find_package(Git)
  execute_process(COMMAND ${GIT_EXECUTABLE} log --pretty=format:%h -n 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE isGit
    OUTPUT_VARIABLE gitHash)
  if(${isGit} EQUAL 0)
    message("  Git commit:             ${gitHash}")
  endif()
endif()
message ("  CMAKE_BUILD_TYPE        ${CMAKE_BUILD_TYPE}")
message ("  CMAKE_C_COMPILER        ${CMAKE_C_COMPILER}")
message ("  CMAKE_CXX_COMPILER      ${CMAKE_CXX_COMPILER}")
message ("  CMAKE_CXX_STANDARD      ${CMAKE_CXX_STANDARD}")

message ("  GISMO_COEFF_TYPE        ${GISMO_COEFF_TYPE}")
message ("  GISMO_INDEX_TYPE        ${GISMO_INDEX_TYPE}")

## #################################################################
## Options list: Standard options
## #################################################################
  
option(GISMO_OPTIONAL_HEAD "Updates submodules to HEAD (online mode)" false  )
if  (${GISMO_OPTIONAL_HEAD})
message ("  GISMO_OPTIONAL_HEAD   ${GISMO_OPTIONAL_HEAD}")
endif()

option(GISMO_PLUGIN_AXL           "Build Axl Plugin"         false  )
if  (${GISMO_PLUGIN_AXL})
message ("  GISMO_PLUGIN_AXL         ${GISMO_PLUGIN_AXL}")
endif()

option(GISMO_BUILD_EXAMPLES      "Build examples"            true   )
if  (${GISMO_BUILD_EXAMPLES})
message ("  GISMO_BUILD_EXAMPLES    ${GISMO_BUILD_EXAMPLES}")
endif()

option(GISMO_BUILD_LIB           "Build shared library"      true   )
message ("  GISMO_BUILD_LIB         ${GISMO_BUILD_LIB}")

option(GISMO_BUILD_MEX           "Build Mex files"           false  )
if  (${GISMO_BUILD_MEX})
message ("  GISMO_BUILD_MEX         ${GISMO_BUILD_MEX}")
endif()

option(GISMO_BUILD_PCH           "Build precompiled headers" false  )
if  (${GISMO_BUILD_PCH})
message ("  GISMO_BUILD_PCH         ${GISMO_BUILD_PCH}")
endif()

option(GISMO_WITH_PYBIND11      "Build Python module using pybind11" false  )
if  (${GISMO_WITH_PYBIND11})
message ("  GISMO_WITH_PYBIND11    ${GISMO_WITH_PYBIND11}")
endif()

option(GISMO_PLUGIN_PVIEW         "Build Paraview Plugin"     false  )
if  (${GISMO_PLUGIN_PVIEW})
message ("  GISMO_PLUGIN_PVIEW      ${GISMO_BUILD_VIEW}")
endif()

option(GISMO_PLUGIN_RHINOS  "Build Rhino THB Plugins"   false  )
if  (${GISMO_PLUGIN_RHINOS})
message ("  GISMO_PLUGIN_RHINOS ${GISMO_PLUGIN_RHINOS}")
endif()

option(GISMO_BUILD_UNITTESTS     "Build unittests"           false  )
if  (${GISMO_BUILD_UNITTESTS})
message ("  GISMO_BUILD_UNITTESTS   ${GISMO_BUILD_UNITTESTS}")
endif()

option(GISMO_WITH_XDEBUG           "Extra debug features"      false  )
if  (${GISMO_WITH_XDEBUG})
message ("  GISMO_WITH_XDEBUG      ${GISMO_EXTRA_XDEBUG}")
endif()

option(GISMO_WITH_ADIFF          "With auto-diff"            false  )
if  (${GISMO_WITH_ADIFF})
message ("  GISMO_WITH_ADIFF        ${GISMO_WITH_ADIFF}")
endif()

#option(GISMO_WITH_METIS          "With METIS"                false )
#if  (${GISMO_WITH_METIS})
#message ("  GISMO_WITH_METIS        ${GISMO_WITH_METIS}")
#endif()

option(GISMO_WITH_MPI            "With MPI"                  false  )
if  (${GISMO_WITH_MPI})
message ("  GISMO_WITH_MPI          ${GISMO_WITH_MPI}")
endif()

option(GISMO_WITH_OPENMP         "With OpenMP"               false  )
if  (${GISMO_WITH_OPENMP})
message ("  GISMO_WITH_OPENMP       ${GISMO_WITH_OPENMP}")
endif()

option(GISMO_WITH_PARDISO        "With PARDISO"              false  )
if  (${GISMO_WITH_PARDISO})
message ("  GISMO_WITH_PARDISO      ${GISMO_WITH_PARDISO}")
endif()

option(GISMO_WITH_PASTIX         "With PastiX"               false  )
if  (${GISMO_WITH_PASTIX})
message ("  GISMO_WITH_PASTIX       ${GISMO_WITH_PASTIX}")
endif()

option(GISMO_WITH_SUPERLU        "With SuperLU"              false  )
if  (${GISMO_WITH_SUPERLU})
message ("  GISMO_WITH_SUPERLU      ${GISMO_WITH_SUPERLU}")
endif()

option(GISMO_WITH_TAUCS          "With Taucs"                false  )
if  (${GISMO_WITH_TAUCS})
message ("  GISMO_WITH_TAUCS        ${GISMO_WITH_TAUCS}")
endif()

option(GISMO_WITH_UMFPACK        "With Umfpack"              false  )
if  (${GISMO_WITH_UMFPACK})
message ("  GISMO_WITH_UMFPACK      ${GISMO_WITH_UMFPACK}")
endif()

## #################################################################
## Options list: Extra options
## #################################################################

option(EIGEN_USE_MKL_ALL         "Eigen use MKL"                 false  )
if (EIGEN_USE_MKL_ALL)
message ("  EIGEN_USE_MKL_ALL       ${EIGEN_USE_MKL_ALL}")
endif()

option(GISMO_BUILD_CPPLOT        "Build cpplot"                  false  )
if (GISMO_BUILD_CPPLOT)
message ("  GISMO_BUILD_CPPLOT      ${GISMO_BUILD_CPPLOT}")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
option(GISMO_BUILD_COVERAGE      "Build with coverage"          false )
if (GISMO_BUILD_COVERAGE)
message ("  GISMO_BUILD_COVERAGE    ${GISMO_BUILD_COVERAGE}")
endif()
endif(CMAKE_COMPILER_IS_GNUCXX)

option(GISMO_BUILD_QT_APP        "Build Qt application"          false  )
if (GISMO_BUILD_QT_APP)
message ("  GISMO_BUILD_QT_APP      ${GISMO_BUILD_QT_APP}")
endif()

option(GISMO_WITH_WARNINGS            "Enable G+Smo related warnings" false  )
if (GISMO_WITH_WARNINGS)
message ("  GISMO_WITH_WARNINGS          ${GISMO_WITH_WARNINGS}")
endif()

option(GISMO_WITH_VTK            "With VTK"                      false  )
if (GISMO_WITH_VTK)
message ("  GISMO_WITH_VTK          ${GISMO_WITH_VTK}")
endif()

if(DEFINED ${isGit} AND ${isGit} EQUAL 0)
  message(STATUS "Type \"${GIT_EXECUTABLE} submodule\" to see the state of submodules")
endif()

#https://www.threadingbuildingblocks.org/documentation
#message ("  GISMO_WITH_ITBB          ${GISMO_WITH_ITBB}")

#message(STATUS "Type cmake -LAH to see all variables")
