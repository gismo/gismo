######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

## #################################################################
## Options list
## #################################################################

#Standard options
option(GISMO_EXTRA_DEBUG         "Extra debug features"   false  )
option(GISMO_BUILD_LIB           "Build shared library"   true   )
#option(GISMO_BUILD_PCH           "Build precompiled headers" false)
option(GISMO_BUILD_EXAMPLES      "Build examples"         true   )
option(GISMO_BUILD_AXL           "Build Axel Plugin"      false  )
option(GISMO_BUILD_PVIEW         "Build Paraview Plugin"  false  )
option(GISMO_BUILD_MEX           "Build Mex files"        false  )
option(GISMO_WITH_OPENMP         "With OpenMP"            false  )
option(GISMO_WITH_PSOLID         "With Parasolid"         false  )
option(GISMO_WITH_MPFR           "With MPFR"              false  )
option(GISMO_WITH_ONURBS         "With OpenNurbs"         false  )
option(GISMO_WITH_IPOPT          "With IpOpt"             false  )
option(GISMO_WITH_SUPERLU        "With SuperLU"           false  )
option(GISMO_WITH_PARDISO        "With PARDISO"           false  )

#Extra options
option(GISMO_BUILD_QT_APP        "Build Qt application"   false  )
option(GISMO_BUILD_CPP11         "Compile using C++11 flags" false)
option(GISMO_WARNINGS            "Enable G+Smo related warnings" false  )
option(GISMO_WITH_VTK            "With VTK"               false  )
option(GISMO_BUILD_CPPLOT        "Build cpplot"           false  )
if(CMAKE_COMPILER_IS_GNUCXX)
option(GISMO_BUILD_COVERAGE      "Build with coverage"    false  )
endif(CMAKE_COMPILER_IS_GNUCXX)

message ("Configuration:")
message ("  Source folder:          ${CMAKE_SOURCE_DIR}")
message ("  CMAKE_BUILD_TYPE        ${CMAKE_BUILD_TYPE}")
message ("  GISMO_COEFF_TYPE        ${GISMO_COEFF_TYPE}")
message ("  GISMO_EXTRA_DEBUG       ${GISMO_EXTRA_DEBUG}")
message ("  GISMO_BUILD_LIB         ${GISMO_BUILD_LIB}")
message ("  GISMO_BUILD_EXAMPLES    ${GISMO_BUILD_EXAMPLES}")
message ("  GISMO_BUILD_AXL         ${GISMO_BUILD_AXL}")
#message ("  GISMO_BUILD_PVIEW       ${GISMO_BUILD_PVIEW}")
#message ("  GISMO_BUILD_MEX         ${GISMO_BUILD_MEX}")
#message ("  GISMO_WITH_OPENMP       ${GISMO_WITH_OPENMP}")
message ("  GISMO_WITH_PSOLID       ${GISMO_WITH_PSOLID}")
#message ("  GISMO_WITH_MPFR         ${GISMO_WITH_MPFR}")
message ("  GISMO_WITH_ONURBS       ${GISMO_WITH_ONURBS}")
#message ("  GISMO_WITH_IPOPT        ${GISMO_WITH_IPOPT}")
message ("  GISMO_WITH_SUPERLU      ${GISMO_WITH_SUPERLU}")

#https://www.threadingbuildingblocks.org/documentation
#message ("  GISMO_WITH_ITBB          ${GISMO_WITH_ITBB}")
