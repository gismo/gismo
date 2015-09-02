
find_library(PARDISO_LIBRARY NAMES pardiso500-GNU481-X86-64 pardiso500-GNU472-X86-64 pardiso500-INTEL1301-X86-64 pardiso500-MACOS-X86-64 libpardiso500-WIN-X86-64
             HINTS ${CMAKE_BINARY_DIR}/lib ${Pardiso_DIR} ${Pardiso_DIR}/lib)

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
