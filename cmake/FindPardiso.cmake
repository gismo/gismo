

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   find_library(PARDISO_LIBRARY NAMES pardiso500-INTEL1301-X86-64 pardiso412-INTEL120-X86-64 pardiso411-INTEL101-X86-64 HINTS ${CMAKE_BINARY_DIR}/lib ${Pardiso_DIR} ${Pardiso_DIR}/lib)
else()
   find_library(PARDISO_LIBRARY NAMES pardiso500-GNU481-X86-64 pardiso500-GNU472-X86-64 pardiso412-GNU450-X86-64 pardiso412-GNU430-X86-64 pardiso411-GNU443-X86-64 pardiso500-MACOS-X86-64 libpardiso500-WIN-X86-64 libpardiso412-WIN-X86-64 libpardiso412-WIN-X86
             HINTS ${CMAKE_BINARY_DIR}/lib ${Pardiso_DIR} ${Pardiso_DIR}/lib)
endif()

if(PARDISO_LIBRARY)
	get_filename_component(PARDISO_LIBRARY_DIR ${PARDISO_LIBRARY} PATH)
   set(PARDISO_FOUND TRUE)
endif()

if(PARDISO_FOUND)
      #message("PARDISO_LIBRARY    : ${PARDISO_LIBRARY}")
   if(NOT PARDISO_FIND_QUIETLY)
      MESSAGE(STATUS "Found Pardiso: ${PARDISO_LIBRARY}")
   endif()
elseif(NOT PARDISO_FOUND)
   set(Pardiso_DIR "Pardiso_DIR-NOTFOUND " CACHE PATH "The path to the PARDISO library")
   if(PARDISO_FIND_REQUIRED)
      message("Could not find PARDISO library.")
      message(FATAL_ERROR "Set the variable Pardiso_DIR and try again.")
   endif()
endif()
