######################################################################
## CMakeLists.txt --- gsFupBasis
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

project(gsFupBasis LANGUAGES CXX Fortran)

enable_language(Fortran)

# Apply same configuration as G+Smo
include(gsConfig)

##
## gsFupBasis extension
##
include_directories(SYSTEM ${OCC_INCLUDE_DIR})
set(FUP_FORT_SOURCES ModuleFup.f90)

# Compile the extension as part of the G+Smo library
add_library(${PROJECT_NAME} OBJECT ${FUP_FORT_SOURCES})

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE Fortran
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions" )

# Add gsFupBasis extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")

#find_library(GFORTRAN_LIBRARY NAMES gfortran)

# Link G+Smo to FupBasis library (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} gfortran
    CACHE INTERNAL "Gismo extra linker objects")

# Install gsFupBasis extension header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsFupBasis/
        FILES_MATCHING PATTERN "*.h")
