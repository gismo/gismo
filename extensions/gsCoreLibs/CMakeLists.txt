######################################################################
### CMakeLists.txt --- gsCoreLibs
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris, Matthias Moller
######################################################################

# Core libraries extension
project(gsCoreLibsExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

##
## BLAS/Lapack library
##
if(GISMO_WITH_SUPERLU OR GISMO_WITH_TAUCS OR GISMO_WITH_UMFPACK)
  find_package(BLAS QUIET)
  find_package(Lapack QUIET)

  if(NOT BLAS_FOUND OR NOT Lapack_FOUND)
    # Set OpenBLAS version
    set(OPENBLAS_VER "0.3.13")

    # Intel
    if(${TARGET_ARCHITECTURE} STREQUAL "atom")
      set(TARGET "ATOM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "bonnell")
      set(TARGET "ATOM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "broadwell")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "cannonlake")
      set(TARGET "SKYLAKEX")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "cascadelake")
      set(TARGET "SKYLAKEX")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "cooperlake")
      set(TARGET "COOPERLAKE")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "core")
      set(TARGET "CORE2")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "dunnigton")
      set(TARGET "DUNNIGTON")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "goldmont")
      set(TARGET "NEHALEM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "goldmont-plus")
      set(TARGET "NEHALEM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "haswell")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "icelake")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "icelake-avx512")
      set(TARGET "SKYLAKEX")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "ivybridge")
      set(TARGET "SANDYBRIDGE")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "kabylake")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "knm")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "knl")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "merom")
      set(TARGET "CORE2")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "nano")
      set(TARGET "NANO")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "nehalem")
      set(TARGET "NEHALEM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "penryn")
      set(TARGET "PENRYN")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "prescott")
      set(TARGET "PRESCOTT")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "sandybridge")
      set(TARGET "SANDYBRIDGE")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "silvermont")
      set(TARGET "NEHALEM")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "skylake")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "skylake-avx512")
      set(TARGET "SKYLAKEX")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "tigerlake")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "tremont")
      set(TARGET "HASWELL")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "westmere")
      set(TARGET "NEHALEM")

    # AMD
    elseif(${TARGET_ARCHITECTURE} STREQUAL "barcelona")
      set(TARGET "BARCELONA")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "bobcat" OR
        ${TARGET_ARCHITECTURE} STREQUAL "amd14h")
      set(TARGET "BOBCAT")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "bulldozer")
      set(TARGET "BULLDOZER")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "excavator")
      set(TARGET "EXCAVATOR")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "k8")
      set(TARGET "OPTERON")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "k8-sse3")
      set(TARGET "OPTERON_SSE3")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "piledriver")
      set(TARGET "PILEDRIVER")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "steamroller")
      set(TARGET "STEAMROLLER")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "zen")
      set(TARGET "ZEN")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "zen2")
      set(TARGET "ZEN")
    elseif(${TARGET_ARCHITECTURE} STREQUAL "zen3")
      set(TARGET "ZEN")

    else()
      set(TARGET "CORE2")
    endif()

    # Download OpenBLAS sources at configure time
    include(gsFetch)
    gismo_fetch_directory(openblas
        URL            https://github.com/xianyi/OpenBLAS/releases/download/v${OPENBLAS_VER}/OpenBLAS-${OPENBLAS_VER}.tar.gz
        DESTINATION    external
        )

    # Build OpenBLAS library at compile time
    add_subdirectory(${gismo_externals}/openblas ${gismo_BINARY_DIR}/external/openblas)

    # Set OpenBLAS library
    set(BLAS_LIBRARIES ${gismo_BINARY_DIR}/external/openblas/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    set(Lapack_LIBRARIES ${gismo_BINARY_DIR}/external/openblas/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    
    find_package(Threads QUIET REQUIRED)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

  endif(NOT BLAS_FOUND OR NOT Lapack_FOUND)
endif(GISMO_WITH_SUPERLU OR GISMO_WITH_TAUCS OR GISMO_WITH_UMFPACK)



##
## SuperLU library
##
if(GISMO_WITH_SUPERLU)
  # Look for pre-installed SuperLU library
  find_package(SuperLU QUIET)
  
  if(NOT SuperLU_FOUND)
    # Set SuperLU version
    set(SUPERLU_VER "4.3")
    
    # Download SuperLU sources at configure time
    include(gsFetch)
    gismo_fetch_directory(superlu
        URL            https://portal.nersc.gov/project/sparse/superlu/superlu_${SUPERLU_VER}.tar.gz
        DESTINATION    external
        )

    # Generate include makefile
    enable_language(Fortran)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/superlu_make.inc ${CMAKE_CURRENT_BINARY_DIR}/superlu_make.inc)

    # Build SuperLU library at compile time
    include(ExternalProject)
    ExternalProject_Add(superlu
        SOURCE_DIR           ${CMAKE_CURRENT_BINARY_DIR}/superlu

        # We use the download step to make a local copy of the source directory and perform in-source build
        DOWNLOAD_COMMAND     ${CMAKE_COMMAND} -E copy_directory ${gismo_externals}/superlu ${CMAKE_CURRENT_BINARY_DIR}/superlu
        
        CONFIGURE_COMMAND    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/superlu_make.inc ${CMAKE_CURRENT_BINARY_DIR}/superlu/make.inc

        BUILD_COMMAND        ${CMAKE_MAKE_PROGRAM} superlulib
        BUILD_IN_SOURCE      1

        UPDATE_COMMAND       ""
        INSTALL_COMMAND      ""
        )
    
    # Set SuperLU library and include directories
    set(SUPERLU_INCLUDES ${gismo_externals}/superlu/SRC CACHE INTERNAL "")
    set(SUPERLU_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/superlu/lib/${CMAKE_STATIC_LIBRARY_PREFIX}superlu_${SUPERLU_VER}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")    
    include_directories(${SUPERLU_INCLUDES})
    
    # Install SuperLU header files
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/superlu/SRC
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h")
  endif(NOT SuperLU_FOUND)
endif(GISMO_WITH_SUPERLU)



##
## TAUCS library
##
if(GISMO_WITH_TAUCS)
  # Look for pre-installed aucs library
  find_package(Taucs QUIET)
  
  if(NOT Taucs_FOUND)
    # Set Taucs version
    set(TAUCS_VER "2.2")  

    # Download Taucs sources at configure time
    include(gsFetch)
    gismo_fetch_directory(taucs
        URL            https://www.tau.ac.il/~stoledo/taucs/${TAUCS_VER}/taucs_full.zip
        DESTINATION    external
        )

    # Build Taucs library at compile time
    include(ExternalProject)
    ExternalProject_Add(taucs1
        SOURCE_DIR           ${CMAKE_CURRENT_BINARY_DIR}/taucs

        # We use the download step to make a local copy of the source directory and perform in-source build
        DOWNLOAD_COMMAND     ${CMAKE_COMMAND} -E copy_directory ${gismo_externals}/taucs ${CMAKE_CURRENT_BINARY_DIR}/taucs
        
        CONFIGURE_COMMAND    ""

        BUILD_COMMAND        ""
        BUILD_IN_SOURCE      1

        UPDATE_COMMAND       ""
        INSTALL_COMMAND      ""
        )

    # Set Taucs library and include directories
    set(TAUCS_INCLUDES ${gismo_externals}/superlu/SRC CACHE INTERNAL "")
    set(TAUCS_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/superlu/lib/${CMAKE_STATIC_LIBRARY_PREFIX}superlu_${SUPERLU_VER}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")    
    include_directories(${TAUCS_INCLUDES})

  endif(NOT Taucs_FOUND)
endif(GISMO_WITH_TAUCS)



##
## Umfpack library
##
if(GISMO_WITH_UMFPACK)
  # Look for pre-installed Umfpack library
  find_package(Umfpack QUIET)
  
  if(NOT Umfpack_FOUND)
    # Set Umfpack version
    set(UMFPACK_VER "4.5.6")
    
    # Download Umfpack sources at configure time
    include(gsFetch)
    gismo_fetch_directory(umfpack
        URL            https://people.engr.tamu.edu/davis/SuiteSparse/SuiteSparse-${UMFPACK_VER}.tar.gz
        DESTINATION    external
        )

    # Build Umfpack library at compile time
    include(ExternalProject)
    ExternalProject_Add(umfpack
        SOURCE_DIR           ${CMAKE_CURRENT_BINARY_DIR}/umfpack

        # We use the download step to make a local copy of the source directory and perform in-source build
        DOWNLOAD_COMMAND     ${CMAKE_COMMAND} -E copy_directory ${gismo_externals}/umfpack ${CMAKE_CURRENT_BINARY_DIR}/umfpack

        BUILD_COMMAND        ${CMAKE_MAKE_PROGRAM} static
        BUILD_IN_SOURCE      1

        CONFIGURE_COMMAND    ""
        UPDATE_COMMAND       ""
        INSTALL_COMMAND      ""
        )
    
    # Set Umfpack libraries
    set(UMFPACK_LIBRARIES 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/AMD/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}amd${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/BTF/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}btf${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/CAMD/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}camd${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/CCOLAMD/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}ccolamd${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/CHOLMOD/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}cholmod${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/COLAMD/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}colamd${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/KLU/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}klu${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/LDL/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}ldl${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/RBio/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}rbio${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/SPQR/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}spqr${CMAKE_STATIC_LIBRARY_SUFFIX} 
      ${CMAKE_CURRENT_BINARY_DIR}/umfpack/UMFPACK/Lib/${CMAKE_STATIC_LIBRARY_PREFIX}umfpack${CMAKE_STATIC_LIBRARY_SUFFIX} 
      CACHE INTERNAL "")
 
    # Set Umfpack include directories
    set(UMFPACK_INCLUDES
      ${gismo_externals}/umfpack/AMD/Include
      ${gismo_externals}/umfpack/BTF/Include
      ${gismo_externals}/umfpack/AMD/Include
      ${gismo_externals}/umfpack/CAMD/Include
      ${gismo_externals}/umfpack/AMD/Include
      ${gismo_externals}/umfpack/CCOLAMD/Include
      ${gismo_externals}/umfpack/CHOLMOD/Include
      ${gismo_externals}/umfpack/COLAMD/Include
      ${gismo_externals}/umfpack/KLU/Include
      ${gismo_externals}/umfpack/LDL/Include
      ${gismo_externals}/umfpack/RBio/Include
      ${gismo_externals}/umfpack/SPQR/Include
      ${gismo_externals}/umfpack/UMFPACK/Include
      CACHE INTERNAL "")
    include_directories(${UMFPACK_INCLUDES})

    # Install Umfpack header files
    install(DIRECTORY ${UMFPACK_INCLUDES}
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h")
  endif(NOT Umfpack_FOUND)
endif(GISMO_WITH_UMFPACK)



##
## gsCoreLibs extension
##

# Compile gsCoreLibs extension as part of the G+Smo library
add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
)

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

# Add gsCoreLibs extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "gismo extensions to be included")

# Add Umfpack include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${BLAS_INCLUDES} ${SUPERLU_INCLUDES} ${TAUCS_INCLUDES} ${UMFPACK_INCLUDES}
  CACHE INTERNAL "gismo include directories")

# Link G+Smo to Umfpack libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${BLAS_LIBRARIES} ${SUPERLU_LIBRARIES} ${TAUCS_LIBRARIES} ${UMFPACK_LIBRARIES}
    CACHE INTERNAL "gismo extra linker objects")

# Install gsCoreLibs header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsCoreLibs/
        FILES_MATCHING PATTERN "*.h")

##
## Add dependencies
##
if(GISMO_WITH_SUPERLU OR GISMO_WITH_TAUCS OR GISMO_WITH_UMFPACK)
  if(NOT BLAS_FOUND OR NOT Lapack_FOUND)
    # Make gsCoreLibs extension dependent on OpenBLAS library
    add_dependencies(${PROJECT_NAME} openblas)
  endif()
endif()

if(GISMO_WITH_SUPERLU AND NOT SuperLU_FOUND)
  # Make gsCoreLibs extension dependent on SuperLU library
  add_dependencies(${PROJECT_NAME} superlu)
  if(NOT BLAS_FOUND OR NOT Lapack_FOUND)
    add_dependencies(superlu openblas)
  endif()
endif()

if(GISMO_WITH_TAUCS AND NOT TAUCS_FOUND)
  # Make gsCoreLibs extension dependent on Taucs library
  add_dependencies(${PROJECT_NAME} taucs)
endif()

if(GISMO_WITH_UMFPACK AND NOT Umfpack_FOUND)
  # Make gsCoreLibs extension dependent on Umfpack library
  add_dependencies(${PROJECT_NAME} umfpack)
  if(NOT BLAS_FOUND OR NOT Lapack_FOUND)
    add_dependencies(superlu openblas)
  endif()
endif()
