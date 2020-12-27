### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Torch extension
project(gsTorchExtension)

# Apply same configuration as G+Smo
include(gsConfig)

if(CMAKE_CXX_STANDARD LESS 14)
  message(FATAL_ERROR "gsTorchExtension requires C++14 or better.")
endif()

macro(pytorch_compile_from_source)
  message(FATAL_ERROR "No binary distribution of Torch available for your operating  system")
endmacro()

find_package(Torch 1.7.0 QUIET)

if ( NOT TORCH_FOUND )

  if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )

    if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
      find_package(CUDA QUIET)
      if(CUDA_FOUND)
        if ("x${CUDA_VERSION}" STREQUAL "x11.0")
          if (DEFINED _GLIBCXX_USE_CXX11_ABI AND NOT ${_GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(Torch-linux-cuda110
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-shared-with-deps-1.7.0%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cuda110 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-linux-cxx11-cuda110
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-cxx11-abi-shared-with-deps-1.7.0%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cxx11-cuda110 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          if (DEFINED _GLIBCXX_USE_CXX11_ABI AND NOT ${_GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(Torch-linux-cuda102
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-shared-with-deps-1.7.0.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cuda102 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-linux-cxx11-cuda102
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-cxx11-abi-shared-with-deps-1.7.0.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cxx11-cuda102 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.1")
          if (DEFINED _GLIBCXX_USE_CXX11_ABI AND NOT ${_GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(Torch-linux-cuda101
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-shared-with-deps-1.7.0%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cuda101 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-linux-cxx11-cuda101
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-cxx11-abi-shared-with-deps-1.7.0%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cxx11-cuda101 CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x9.2")
          if (DEFINED _GLIBCXX_USE_CXX11_ABI AND NOT ${_GLIBCXX_USE_CXX11_ABI})
            gismo_fetch_directory(Torch-linux-cuda92
              URL            https://download.pytorch.org/libtorch/cu92/libtorch-shared-with-deps-1.7.0%2Bcu92.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cuda92 CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-linux-cxx11
              URL            https://download.pytorch.org/libtorch/cu92/libtorch-cxx11-abi-shared-with-deps-1.7.0%2Bcu92.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cxx11-cuda92 CACHE INTERNAL "")
          endif()
        else()
          pytorch_compile_trom_source()
        endif()
      else(CUDA_FOUND)
        if (DEFINED _GLIBCXX_USE_CXX11_ABI AND NOT ${_GLIBCXX_USE_CXX11_ABI}) 
          gismo_fetch_directory(Torch-linux
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.7.0%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux CACHE INTERNAL "")
        else()
          gismo_fetch_directory(Torch-linux-cxx11
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.7.0%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-linux-cxx11 CACHE INTERNAL "")
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
        gismo_fetch_directory(Torch-macos
          URL            https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.7.0.zip
          DESTINATION    external
          )
        set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-macos CACHE INTERNAL "")
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
            gismo_fetch_directory(Torch-win-cuda110-release
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-win-shared-with-deps-1.7.0%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda110-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-win-cuda110-debug
              URL            https://download.pytorch.org/libtorch/cu110/libtorch-win-shared-with-deps-debug-1.7.0%2Bcu110.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda110-debug CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.2")
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(Torch-win-cuda102-release
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-1.7.0%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda102-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-win-cuda102-debug
              URL            https://download.pytorch.org/libtorch/cu102/libtorch-win-shared-with-deps-debug-1.7.0%2Bcu102.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda102-debug CACHE INTERNAL "")
          endif()
        elseif ("x${CUDA_VERSION}" STREQUAL "x10.1")
          if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
            gismo_fetch_directory(Torch-win-cuda101-release
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-win-shared-with-deps-1.7.0%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda101-release CACHE INTERNAL "")
          else()
            gismo_fetch_directory(Torch-win-cuda101-debug
              URL            https://download.pytorch.org/libtorch/cu101/libtorch-win-shared-with-deps-debug-1.7.0%2Bcu101.zip
              DESTINATION    external
              )
            set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-cuda101-debug CACHE INTERNAL "")
          endif()
        else()
          pytorch_compile_from_source()
        endif()
      else(CUDA_FOUND)
        if ("x${CMAKE_BUILD_TYPE}" STREQUAL "xRelease")
          gismo_fetch_directory(Torch-win-release
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-1.7.0%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-release CACHE INTERNAL "")
        else()
          gismo_fetch_directory(Torch-win-debug
            URL            https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-debug-1.7.0%2Bcpu.zip
            DESTINATION    external
            )
          set(TORCH_INSTALL_PREFIX ${gismo_SOURCE_DIR}/external/Torch-win-debug CACHE INTERNAL "")
        endif()
      endif(CUDA_FOUND)
    else()
      pytorch_compile_from_source()
    endif()

  else()
    pytorch_compile_from_source()
  endif()
endif ( NOT TORCH_FOUND )

find_library(C10_LIBRARY        NAMES c10        HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_LIBRARY      NAMES torch      HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CPU_LIBRARY  NAMES torch_cpu  HINTS "${TORCH_INSTALL_PREFIX}/lib")
find_library(TORCH_CUDA_LIBRARY NAMES torch_cuda HINTS "${TORCH_INSTALL_PREFIX}/lib")

set(TORCH_LIBRARIES ${TORCH_LIBRARY} ${TORCH_CPU_LIBRARY} ${C10_LIBRARY} CACHE INTERNAL "")
set(TORCH_INCLUDE_DIRS ${TORCH_INSTALL_PREFIX}/include ${TORCH_INSTALL_PREFIX}/include/torch/csrc/api/include CACHE INTERNAL "")

include_directories(${TORCH_INCLUDE_DIRS})

set(TORCH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR} gsTorch_.cpp)
add_library(${PROJECT_NAME} OBJECT ${TORCH_SOURCES})

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"  )

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "Gismo extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${TORCH_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${TORCH_INCLUDE_DIRS}
  DESTINATION include/gismo/gsTorch/
  FILES_MATCHING PATTERN "*.h")
