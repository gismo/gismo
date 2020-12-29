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
set(TORCH_VER 1.7.1)

# Apply same configuration as G+Smo
include(gsConfig)

# Check that C++14 is enabled
if(CMAKE_CXX_STANDARD LESS 14)
  message(FATAL_ERROR "gsTorchExtension requires C++14 or better.")
endif()

# Compile LibTorch from source
macro(pytorch_compile_from_source)
  message(FATAL_ERROR "No binary distribution of Torch available for your operating  system")
endmacro()

# Look for pre-installed LibTorch libraries
find_package(Torch ${TORCH_VER} QUIET)

if ( NOT TORCH_FOUND )

  if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )

    # Compile test program to check whether pre-cxx11 or cxx11 ABI is used
    try_compile(GLIBCXX_USE_CXX11_ABI ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test_linux_cxx11_abi.cxx)

    # When we build libtorch with the old GCC ABI, dependent libraries must too.
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if (${GLIBCXX_USE_CXX11_ABI})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
      else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=1")
      endif()
    endif()
    
    if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      find_package(CUDA QUIET)
      if(CUDA_FOUND)
        # CUDA 11.0
        if ("x${CUDA_VERSION}" STREQUAL "x11.0")
          if (${GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(torch-linux-cxx11-cuda110
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cxx11-cuda110 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-linux-cuda110
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-shared-with-deps-${TORCH_VER}%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda110 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          # CUDA 10.2
          if (${GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(torch-linux-cxx11-cuda102
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cxx11-cuda102 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-linux-cuda102
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-shared-with-deps-${TORCH_VER}.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda102 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.1")
          # CUDA 10.1
          if (${GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(torch-linux-cxx11-cuda101
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cxx11-cuda101 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-linux-cuda101
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-shared-with-deps-${TORCH_VER}%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda101 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x9.2")
          # CUDA 9.2
          if (${GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(torch-linux-cxx11-cuda92
              URL            https://download.pytorch.org/libtorch/cu92/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}%2Bcu92.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cxx11-cuda92 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-linux-cuda92
              URL            https://download.pytorch.org/libtorch/cu92/libtorch-shared-with-deps-${TORCH_VER}%2Bcu92.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cuda92 CACHE INTERNAL "")
          endif()
        else()
          pytorch_compile_trom_source()
        endif()
      else(CUDA_FOUND)
        # No CUDA
        if (${GLIBCXX_USE_CXX11_ABI})
          gismo_fetch_directory(torch-linux-cxx11
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux-cxx11 CACHE INTERNAL "")
        else()
          gismo_fetch_directory(torch-linux
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-linux CACHE INTERNAL "")
        endif()
      endif(CUDA_FOUND)
    else()
      pytorch_compile_from_source()
    endif()

  elseif ( CMAKE_SYSTEM_NAME MATCHES "Darwin" )

    if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      find_package(CUDA QUIET)
      if(CUDA_FOUND)
        pytorch_compile_from_source()
      else(CUDA_FOUND)
        gismo_fetch_directory(torch-macos
          URL            https://download.pytorch.org/libtorch/cpu/libtorch-macos-${TORCH_VER}.zip
          DESTINATION    external
          )
        set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-macos CACHE INTERNAL "")
      endif(CUDA_FOUND)
    else()
      pytorch_compile_from_source()
    endif()

  elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )

    if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      find_package(CUDA QUIET)
      if(CUDA_FOUND)
        if ("x${CUDA_VERSION}" STREQUAL "x11.0")
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(torch-win-cuda110-release
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda110-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-win-cuda110-debug
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda110-debug CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(torch-win-cuda102-release
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda102-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-win-cuda102-debug
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda102-debug CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.1")
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(torch-win-cuda101-release
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda101-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(torch-win-cuda101-debug
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-cuda101-debug CACHE INTERNAL "")
          endif()
        else()
          pytorch_compile_from_source()
        endif()
      else(CUDA_FOUND)
        if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
          gismo_fetch_directory(torch-win-release
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-release CACHE INTERNAL "")
        else()
          gismo_fetch_directory(torch-win-debug
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-debug-${TORCH_VER}%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/torch-win-debug CACHE INTERNAL "")
        endif()
      endif(CUDA_FOUND)
    else()
      pytorch_compile_from_source()
    endif()

  else()
    pytorch_compile_from_source()
  endif()
endif ( NOT TORCH_FOUND )

# Find LibTorch libraries
find_library(C10_LIBRARY        NAMES c10        HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_LIBRARY      NAMES torch      HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CPU_LIBRARY  NAMES torch_cpu  HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CUDA_LIBRARY NAMES torch_cuda HINTS "${TORCH_INSTALL_PREFIX}/lib")

# Set LibTorch libraries
set(TORCH_LIBRARIES ${TORCH_LIBRARY} ${C10_LIBRARY} ${TORCH_CPU_LIBRARY} ${TORCH_CUDA_LIBRARY} CACHE INTERNAL "")

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
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

# Add gsTorch extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")

# Add LibTorch include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${TORCH_INCLUDE_DIRS}
  CACHE INTERNAL "Gismo include directories")

# Link G+Smo to LibTorch libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${TORCH_LIBRARIES}
    CACHE INTERNAL "Gismo extra linker objects")

# Install gsTorch header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsTorch/
        FILES_MATCHING PATTERN "*.h")
