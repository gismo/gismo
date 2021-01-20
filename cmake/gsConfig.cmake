######################################################################
## gsConfig.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

include(CheckCXXCompilerFlag)

#find_package(Metis REQUIRED)

#Remove NDEBUG from RelWithDebInfo builds
string(REPLACE "-DNDEBUG" "" replacementFlags "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${replacementFlags} CACHE INTERNAL "" FORCE)
string(REPLACE "-DNDEBUG" "" replacementFlags "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO ${replacementFlags} CACHE INTERNAL "" FORCE)

if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  #fixme: enable for Darwin/clang (probably no export explicit template instantiations)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
  set(CMAKE_C_VISIBILITY_PRESET   hidden)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN 1 )
endif()

# Set a default coefficient numeric types if not specified
if(NOT GISMO_COEFF_TYPE)
  set (GISMO_COEFF_TYPE "double" CACHE STRING
   "Coefficient type(float, double, long double, mpfr::mpreal, mpq_class, posit_2_0, posit_3_0, posit_3_1, posit_4_0, posit_8_0, posit_8_1, posit_16_1, posit_32_2, posit_64_3, posit_128_4, posit_256_5)" FORCE)
elseif(${GISMO_COEFF_TYPE} STREQUAL "mpfr::mpreal")
  set(GISMO_WITH_MPFR ON CACHE BOOL "Use MPFR" FORCE)
elseif(${GISMO_COEFF_TYPE} STREQUAL "mpq_class")
  set(GISMO_WITH_GMP ON CACHE BOOL "Use GMP" FORCE)
elseif(${GISMO_COEFF_TYPE} STREQUAL "posit_2_0"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_3_0"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_3_1"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_4_0"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_8_0"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_8_1"   OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_16_1"  OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_32_2"  OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_64_3"  OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_128_4" OR
       ${GISMO_COEFF_TYPE} STREQUAL "posit_256_5")
  set(GISMO_WITH_UNUM ON CACHE BOOL "Use UNUM" FORCE)
endif()
set_property(CACHE GISMO_COEFF_TYPE PROPERTY STRINGS
"float" "double" "long double" "mpfr::mpreal" "mpq_class" "posit_2_0" "posit_3_0" "posit_3_1" "posit_4_0" "posit_8_0" "posit_8_1" "posit_16_1" "posit_32_2" "posit_64_3" "posit_128_4" "posit_256_5")

if(NOT GISMO_INDEX_TYPE)
   set (GISMO_INDEX_TYPE "int" CACHE STRING
   #math(EXPR BITSZ_VOID_P "8*${CMAKE_SIZEOF_VOID_P}")
   #set (GISMO_INDEX_TYPE "int${BITSZ_VOID_P}_t" CACHE STRING
   "Index type(int, int32_t, int64_t, long, long long)" FORCE)
   set_property(CACHE GISMO_INDEX_TYPE PROPERTY STRINGS
   "int" "int32_t" "int64_t" "long" "long long" )
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
   "Type of build (None Debug Release RelWithDebInfo MinSizeRel)" FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
     "RelWithDebInfo" "MinSizeRel")
endif()

set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)
foreach(config ${CMAKE_CONFIGURATION_TYPES}) # For Visual studio
    # overrides Debug/Release subfolders
    string(TOUPPER ${config} CONFIG)
    set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} lib)
    set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} bin)
    set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} lib)
endforeach()

#Configure Valgrind
#find_program( MEMORYCHECK_COMMAND valgrind )
#--gen-suppressions=all --trace-children=yes --track-origins=yes
#set( MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-reachable=yes" CACHE INTERNAL "")
set( MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --leak-check=yes -q" CACHE INTERNAL "") #note: empty defaults to "-q --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50
set( MEMORYCHECK_SUPPRESSIONS_FILE "${gismo_SOURCE_DIR}/cmake/valgrind_supp.txt" CACHE INTERNAL "")

set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS OFF)
include(AddCXXCompileOptions)

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntel")
  # message(STATUS "Using Boost for smart pointers")
  find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# Print compilation statistics (these flags work on GCC compiler only)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-report")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Q")

if(GISMO_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  # see http://www.cmake.org/Wiki/CTest:Coverage
  # and http://cmake.3232098.n2.nabble.com/Running-coverage-analysis-td7145452.html
  include(CodeCoverage)
  APPEND_COVERAGE_COMPILER_FLAGS()
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
  #set(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
endif(GISMO_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")

    #if(${MSVC_VERSION} EQUAL 1800 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.00.31101.0)
    #   message(WARNING "Visual Studio 2013 without Update 4 detected. Update your compiler to avoid G+Smo compilation problems.")
    #endif()

    # Disable checked iterators and irrelevant warnings
    #wd4351: regards old behaviour before MSVC2005
    set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS}  /bigobj /D_SECURE_SCL=0  /wd4351 /MP")
    # See http://msdn.microsoft.com/en-us/library/hh697468.aspx
    #add_definitions(-D_HAS_ITERATOR_DEBUGGING=0)
    #add_definitions(-D_SECURE_SCL=0)
    #add_definitions(-D_ITERATOR_DEBUG_LEVEL=0) #VS2012

    # disable incremental linking for executables (it doesn't help for linking with libraries) -- check
    #STRING(REPLACE "/INCREMENTAL:YES" "/INCREMENTAL:NO" CMAKE_EXE_LINKER_FLAGS_DEBUG ${CMAKE_EXE_LINKER_FLAGS_DEBUG})
    #STRING(REPLACE "/INCREMENTAL:YES" "/INCREMENTAL:NO" CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO})

#    if ( GISMO_BUILD_LIB )
#    # /MD /MDd
#      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")
#    endif()

    if (CMAKE_SIZEOF_VOID_P EQUAL 8) #64bit compiler
       # Note: On 64bit-platforms, /Wp64 flag is present, causing extra warnings
       set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} /wd4244 /wd4267")

    #else() #32bit compiler has CMAKE_SIZEOF_VOID_P EQUAL 4
    endif()

endif()

if(GISMO_EXTRA_DEBUG)
  include(gsDebugExtra)
endif(GISMO_EXTRA_DEBUG)

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  else()
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /W4")
  endif()

elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Note: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431
  # affects -Wno-ignored-attributes in Eigen
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wunused-variable") # -fmax-errors=5
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") #-ftrack-macro-expansion=0 -Wno-ignored-attributes
  endif()
  if ("x${CMAKE_CXX_STANDARD}" STREQUAL "x98"
      AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.4)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++11-compat")
  endif()

  if(GISMO_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual -Wextra")
    #-Wshadow -Wconversion -pedantic -Wunused -Wattributes
  endif()

endif()

if (CMAKE_COMPILER_IS_GNUCXX AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  #-Wl,--no-allow-shlib-undefined
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
  if (NOT MINGW)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,defs")
  endif()
endif()

if (MINGW)
  # fixme: export explicit template instantiations in MinGW ?
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--export-all-symbols")

  # large files can overflow pe/coff sections, so use the pe+ format
  CHECK_CXX_COMPILER_FLAG("-Wa,-mbig-obj" HAS_MBIGOBJ)
  if(NOT HAS_MBIGOBJ)
    message(WARNING "Current compiler does not suppport -Wa,-mbig-obj option.")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffunction-sections -Wl,--gc-sections")
  endif()
elseif(NOT MSVC AND NOT POLICY CMP0063 AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  #fixme: enable for Darwin/clang (probably no export explicit template instantiations)
  check_cxx_compiler_flag(-fvisibility=hidden visibility)
    if (visibility) # for object libraries with cmake less than 3.3
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden")
    endif()
endif()

if (GISMO_WITH_OPENMP)
   find_package(OpenMP REQUIRED)
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (CMAKE_COMPILER_IS_GNUCXX AND NOT GISMO_WITH_OPENMP)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT GISMO_WITH_OPENMP)
   if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 3180") #comma for more warns
   elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qdiag-disable:3180")
   endif()
   #set_property(TARGET mytarget PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
endif()

#CHECK_CXX_SOURCE_COMPILES(
#"template<typename T> class A {}; extern template class A<int>; int main() {}"
#GISMO_HAS_EXTERN_TEMPLATES)

#message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")
#message("CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
#message("CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
#message("CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
#message(STATUS "Using compilation flags: ${CMAKE_CXX_FLAGS}, ${CMAKE_CXX_FLAGS_${TEMP}}")

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  #https://github.com/VcDevel/Vc/blob/master/cmake/OptimizeForArchitecture.cmake
  include( OptimizeForArchitecture )
  OptimizeForArchitecture()
  foreach (flag ${OFA_ARCHITECTURE_FLAGS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  endforeach()
endif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
