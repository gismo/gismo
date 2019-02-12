
include(gsFetch)
gismo_fetch_module(devel)

# Fetching surface_mesh external

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

gismo_fetch_directory(surface_mesh
  GIT_REPOSITORY https://opensource.cit-ec.de/git/surface_mesh
  # todo: tarball etc
)

# Case 1: Use third party CMake configuration
#add_subdirectory(${gismo_SOURCE_DIR}/extensions/surface_mesh ${gismo_BINARY_DIR}/extensions/surface_mesh EXCLUDE_FROM_ALL)


# Case 2: Use our own cmake script
file(GLOB SMSOURCES ${gismo_SOURCE_DIR}/extensions/surface_mesh/src/surface_mesh/*.cpp)
file(GLOB SMHEADERS ${gismo_SOURCE_DIR}/extensions/surface_mesh/src/surface_mesh/*.h)
set(SMESH_INCLUDE_DIRS ${gismo_SOURCE_DIR}/extensions/surface_mesh/src CACHE INTERNAL "surface_mesh include dir")
#.. target_include_directories
add_library(surface_mesh OBJECT EXCLUDE_FROM_ALL ${SMSOURCES} ${SMHEADERS})

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:surface_mesh>
    CACHE INTERNAL "Gismo extensions to be included")

#install(TARGETS surface_mesh DESTINATION lib OPTIONAL)
install( FILES ${HEADERS} DESTINATION include/surface_mesh/ OPTIONAL)

# Case 3: use ExternalProject_Add (to do)
