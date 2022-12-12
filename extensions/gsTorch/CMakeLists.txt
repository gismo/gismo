### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Torch extension
project(gsTorchExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Set LibTorch version
set(TORCH_VER 1.10.1)
set(TORCH_MINVER 1.8.2)

# Apply same configuration as G+Smo
include(gsConfig)

# Check that C++14 is enabled and that compilers are not too old
if (CMAKE_CXX_STANDARD LESS 14 OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"        AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5) OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"      AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4) OR
    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4) )
  message(FATAL_ERROR "gsTorchExtension requires C++14 and GCC 5 or Clang 4 or better.")
endif()

# Compile LibTorch from source
macro(pytorch_compile_from_source)
  message(STATUS "No binary distribution of LibTorch available for your operating system. Building LibTorch from source.")
  include(gsFetch)
  gismo_fetch_directory(torch-src
    URL            https://github.com/pytorch/pytorch/releases/download/v${TORCH_VER}/pytorch-v${TORCH_VER}.tar.gz
    DESTINATION    external
    )
  set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-src CACHE INTERNAL "")  
endmacro()

# Look for pre-installed LibTorch libraries
find_package(Torch ${TORCH_MINVER} QUIET)

if ( NOT TORCH_FOUND )

  include(gsFetch)

  if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    ##
    ## Linux
    ##

    if ("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      ##
      ## x86/x86_64
      ##
      
      find_package(CUDA QUIET)
      if (CUDA_FOUND)

        if ("x${CUDA_VERSION}" STREQUAL "x11.3")
          ##
          ## CUDA 11.3
          ##
          gismo_fetch_directory(torch-linux-cuda113
            URL            https://download.pytorch.org/libtorch/cu113/libtorch-shared-with-deps-${TORCH_VER}%2Bcu113.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda113 CACHE INTERNAL "")

        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          ##
          ## CUDA 10.2
          ##
          gismo_fetch_directory(torch-linux-cuda102
            URL            https://download.pytorch.org/libtorch/cu102/libtorch-shared-with-deps-${TORCH_VER}.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda102 CACHE INTERNAL "")

        else()
          ##
          ## CUDA other version
          ##
          pytorch_compile_trom_source()
        endif()

      else()
        ##
        ## No CUDA
        ##
        gismo_fetch_directory(torch-linux
          URL            https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-${TORCH_VER}%2Bcpu.zip
          DESTINATION    external
          )
        set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux CACHE INTERNAL "")
      endif()

    else()
      ##
      ## Linux non-x86/x86_64
      ##
      pytorch_compile_from_source()
    endif()

  elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    ##
    ## macOS (aka Darwin)
    ##

    if ("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      ##
      ## x86/x86_64
      ##
      
      find_package(CUDA QUIET)
      if (CUDA_FOUND)
        ##
        ## CUDA
        ##
        pytorch_compile_from_source()

      else()
        ##
        ## No CUDA
        ##
        gismo_fetch_directory(torch-macos
          URL            https://download.pytorch.org/libtorch/cpu/libtorch-macos-${TORCH_VER}.zip
          DESTINATION    external
          )
        set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-macos CACHE INTERNAL "")
      endif()

    else()
      ##
      ## non-x86/x86_64
      ##
      pytorch_compile_from_source()
    endif()

  elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
    ##
    ## Windows
    ##

    if ("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      ##
      ## x86/x86_64
      ##
      
      find_package(CUDA QUIET)
      if (CUDA_FOUND)

        if ("x${CUDA_VERSION}" STREQUAL "x11.3")
          ##
          ## CUDA 11.3
          ##
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(torch-win-cuda113-release
              URL            https://download.pytorch.org/libtorch/cu113/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcu113.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda113-release CACHE INTERNAL "")
          elseif ("x${CMAKE_BUILD_TYPE}" STREQUAL "xDebug")
            gismo_fetch_directory(torch-win-cuda113-debug
              URL            https://download.pytorch.org/libtorch/cu113/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcu113.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda113-debug CACHE INTERNAL "")
          else()
            pytorch_compile_from_source()
          endif()
          
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          ##
          ## CUDA 10.2
          ##
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(torch-win-cuda102-release
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda102-release CACHE INTERNAL "")
          elseif ("x${CMAKE_BUILD_TYPE}" STREQUAL "xDebug")
            gismo_fetch_directory(torch-win-cuda102-debug
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda102-debug CACHE INTERNAL "")
          else()
            pytorch_compile_from_source()  
          endif()

        else()
          ##
          ## CUDA other version
          ##
          pytorch_compile_from_source()
        endif()

      else()
        ##
        ## No CUDA
        ##
        if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
          gismo_fetch_directory(torch-win-release
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-release CACHE INTERNAL "")
        elseif ("x${CMAKE_BUILD_TYPE}" STREQUAL "xDebug")
          gismo_fetch_directory(torch-win-debug
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-debug CACHE INTERNAL "")
        else()
          pytorch_compile_from_source()
        endif()
      endif()

    else()
      ##
      ## non-x86/x86_64
      ##
      pytorch_compile_from_source()
    endif()

  else()
    ##
    ## non-standard operating system
    ##
    pytorch_compile_from_source()
  endif()
endif()

# Find LibTorch libraries
find_library(C10_LIBRARY        NAMES c10        HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_LIBRARY      NAMES torch      HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CPU_LIBRARY  NAMES torch_cpu  HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CUDA_LIBRARY NAMES torch_cuda HINTS "${TORCH_INSTALL_PREFIX}/lib")

# Set LibTorch libraries
if (TORCH_CUDA_LIBRARY)
  set(TORCH_LIBRARIES ${TORCH_LIBRARY} ${C10_LIBRARY} ${TORCH_CPU_LIBRARY} ${TORCH_CUDA_LIBRARY} CACHE INTERNAL "")
else()
  set(TORCH_LIBRARIES ${TORCH_LIBRARY} ${C10_LIBRARY} ${TORCH_CPU_LIBRARY} CACHE INTERNAL "")
endif()

# Set LibTorch include directories
set(TORCH_INCLUDE_DIRS ${TORCH_INSTALL_PREFIX}/include ${TORCH_INSTALL_PREFIX}/include/torch/csrc/api/include CACHE INTERNAL "")
include_directories(${TORCH_INCLUDE_DIRS})

##
## gsTorch extension
##

# Compile gsTorch extension as part of the G+Smo library
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
    #START Export all symbols from this extension
    CXX_VISIBILITY_PRESET default
    C_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN 0
    #END Export all symbols from this extension
    FOLDER "G+Smo extensions"
    )

# Add gsTorch extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")

# Add LibTorch include directories to G+Smo standard include directories
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TORCH_INCLUDE_DIRS}
    CACHE INTERNAL "Gismo include directories")

# Link G+Smo to LibTorch libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${TORCH_LIBRARIES}
    CACHE INTERNAL "Gismo extra linker objects")

# Install gsTorch header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsTorch/
        FILES_MATCHING PATTERN "*.h")
      
# Add filedata folder
add_definitions(-DTORCH_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/") 

# Add example files
include_directories(${CODIPACK_INCLUDE_DIR})
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
  add_gismo_executable(${file})
  get_filename_component(tarname ${file} NAME_WE) # name without extension
  set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
  set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
  if( GISMO_WITH_MPI )
    target_include_directories(${tarname} PRIVATE ${MPI_INCLUDE_PATH})
  endif()
  # Install the example executables (optionally)
  install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

# Add unittest files
if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(unittests)
endif(GISMO_BUILD_UNITTESTS)
