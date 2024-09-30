# FIXME: How do I find the version of PARASOLID that I want to use?
# What versions are available?

# NOTE: Parasoid_DIR is understood to be the path to the root of the PARASOLID
# installation library.
set(Parasolid_DIR "PARASOLID_DIR-NOTFOUND " CACHE PATH "The path to the PARASOLID library")

find_path(PARASOLID_INCLUDE_DIR parasolid_kernel.h #parasolid_typedefs.h 
	PATHS ${Parasolid_DIR} ${Parasolid_DIR}/shared_object
)

find_library(PARASOLID_LIBRARY NAMES pskernel
	PATHS ${Parasolid_DIR}/shared_object
)

if(PARASOLID_INCLUDE_DIR AND PARASOLID_LIBRARY)
	get_filename_component(PARASOLID_LIBRARY_DIR ${PARASOLID_LIBRARY} PATH)
   set(Parasolid_FOUND TRUE)
endif()


if(Parasolid_FOUND)
   if(NOT Parasolid_FIND_QUIETLY)
      MESSAGE(STATUS "Found Parasolid: ${PARASOLID_LIBRARY}")
   endif()
elseif(NOT Parasolid_FOUND)
   if(NOT Parasolid_FIND_QUIETLY)
      message("Could not find Parasolid library.")
   endif()
  if(Parasolid_FIND_REQUIRED)
      message(FATAL_ERROR "Set the variable Parasolid_DIR and try again.")
   endif(Parasolid_FIND_REQUIRED)
endif()
