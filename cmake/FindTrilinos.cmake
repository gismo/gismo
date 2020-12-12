######################################################################
## FindTrilinos.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
##
## Required Trilinos packages are:
## - Epetra
## - Teuchos
##
######################################################################

unset(Trilinos_LIBRARIES CACHE)
#unset(AZTECOO_FOUND CACHE)
#unset(AZTECOO_LIBRARIES CACHE)
  
FIND_PATH(AMESOS_INCLUDE_PATH                 Amesos.h                                                 PATHS ${Trilinos_DIR}/include)
FIND_PATH(AZTECOO_INCLUDE_PATH                AztecOO.h                                                PATHS ${Trilinos_DIR}/include)
FIND_PATH(EPETRA_INCLUDE_PATH                 Epetra_Object.h                                          PATHS ${Trilinos_DIR}/include)
FIND_PATH(IFPACK_INCLUDE_PATH                 Ifpack.h                                                 PATHS ${Trilinos_DIR}/include)
FIND_PATH(LOCA_INCLUDE_PATH                   LOCA.H                                                   PATHS ${Trilinos_DIR}/include)
FIND_PATH(ML_INCLUDE_PATH                     MLAPI.h                                                  PATHS ${Trilinos_DIR}/include)
FIND_PATH(NOX_INCLUDE_PATH                    NOX.H                                                    PATHS ${Trilinos_DIR}/include)
FIND_PATH(TEUCHOS_INCLUDE_PATH                Teuchos_Object.hpp                                       PATHS ${Trilinos_DIR}/include)
FIND_PATH(KOMPLEX_INCLUDE_PATH                Komplex_Version.h                                        PATHS ${Trilinos_DIR}/include)

FIND_PATH(LOCA_EPETRA_INCLUDE_PATH            LOCA_Epetra.H                                            PATHS ${Trilinos_DIR}/include)
FIND_PATH(NOX_EPETRA_INCLUDE_PATH             NOX_Epetra.H                                             PATHS ${Trilinos_DIR}/include)
FIND_PATH(EPETRAEXT_INCLUDE_PATH              EpetraExt_Version.h                                      PATHS ${Trilinos_DIR}/include)

FIND_LIBRARY(AMESOS_LIBRARY                   NAMES amesos trilinos_amesos                             PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(AZTECOO_LIBRARY                  NAMES aztecoo trilinos_aztecoo                           PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(EPETRA_LIBRARY                   NAMES epetra trilinos_epetra                             PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(IFPACK_LIBRARY                   NAMES ifpack trilinos_ifpack                             PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(LOCA_LIBRARY                     NAMES loca trilinos_loca                                 PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(ML_LIBRARY                       NAMES ml trilinos_ml                                     PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(NOX_LIBRARY                      NAMES nox trilinos_nox                                   PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(KOKKOS_LIBRARY                   NAMES kokkos trilinos_kokkos kokkoscore                  PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOSKOKKOSCOMM_LIBRARY        NAMES teuchoskokkoscomm trilinos_teuchoskokkoscomm       PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOSKOKKOSCOMPAT_LIBRARY      NAMES teuchoskokkoscompat trilinos_teuchoskokkoscompat   PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOS_LIBRARY                  NAMES teuchos trilinos_teuchos teuchoscore               PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOS_LIBRARY_COMM             NAMES teuchoscomm                                        PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOS_LIBRARY_NUMERICS         NAMES teuchosnumerics                                    PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOS_LIBRARY_PARAMETER_LIST   NAMES teuchosparameterlist                               PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(TEUCHOS_LIBRARY_REMAINDER        NAMES teuchosremainder                                   PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(KOMPLEX_LIBRARY                  NAMES komplex trilinos_komplex                           PATHS ${Trilinos_DIR}/lib)

FIND_LIBRARY(THYRA_LIBRARY                    NAMES thyra thyracore                                    PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(THYRA_EPETRA_LIBRARY             NAMES thyraepetra trilinos_thyraepetra                   PATHS ${Trilinos_DIR}/lib)

FIND_LIBRARY(LOCA_EPETRA_LIBRARY              NAMES locaepetra trilinos_locaepetra                     PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(NOX_EPETRA_LIBRARY               NAMES noxepetra trilinos_noxepetra                       PATHS ${Trilinos_DIR}/lib)
FIND_LIBRARY(EPETRAEXT_LIBRARY                NAMES epetraext trilinos_epetraext                       PATHS ${Trilinos_DIR}/lib)

# Experimental
if(WITH_ZOLTAN)
  FIND_PATH(ZOLTAN_INCLUDE_PATH               zoltan.h                                                 PATHS ${Trilinos_DIR}/include)
  FIND_LIBRARY(ZOLTAN_LIBRARY                 zoltan                                                   PATHS ${Trilinos_DIR}/lib)
endif(WITH_ZOLTAN)

INCLUDE(FindPackageHandleStandardArgs)

IF(EPETRA_INCLUDE_PATH AND EPETRA_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${EPETRA_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${EPETRA_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_EPETRA YES)
ENDIF(EPETRA_INCLUDE_PATH AND EPETRA_LIBRARY)

IF(TEUCHOS_INCLUDE_PATH AND TEUCHOS_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${TEUCHOS_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES}  ${TEUCHOS_LIBRARY_NUMERICS} ${TEUCHOS_LIBRARY_REMAINDER} ${TEUCHOS_LIBRARY_COMM} ${TEUCHOS_LIBRARY} ${TEUCHOS_LIBRARY_PARAMETER_LIST}
   CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_TEUCHOS YES)
ENDIF(TEUCHOS_INCLUDE_PATH AND TEUCHOS_LIBRARY)

IF(KOKKOS_LIBRARY)
  #SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${KOKKOS_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${KOKKOS_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_KOKKOS YES)
ENDIF()

IF(TEUCHOSKOKKOSCOMPAT_LIBRARY)
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${TEUCHOSKOKKOSCOMPAT_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_KOKKOSCOMPAT YES)
ENDIF()

IF(TEUCHOSKOKKOSCOMM_LIBRARY)
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${TEUCHOSKOKKOSCOMM_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_KOKKOSCOMM YES)
ENDIF()

find_package_handle_standard_args(EPETRA DEFAULT_MSG EPETRA_LIBRARY)
find_package_handle_standard_args(TEUCHOS DEFAULT_MSG TEUCHOS_LIBRARY)

IF(AMESOS_INCLUDE_PATH AND AMESOS_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${AMESOS_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${AMESOS_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_AMESOS YES)
  find_package_handle_standard_args(AMESOS DEFAULT_MSG AMESOS_LIBRARY)
ENDIF(AMESOS_INCLUDE_PATH AND AMESOS_LIBRARY)

IF(AZTECOO_INCLUDE_PATH AND AZTECOO_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${AZTECOO_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${AZTECOO_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_AZTECOO YES)
  find_package_handle_standard_args(AZTECOO DEFAULT_MSG AZTECOO_LIBRARY)
ENDIF(AZTECOO_INCLUDE_PATH AND AZTECOO_LIBRARY)

IF(IFPACK_INCLUDE_PATH AND IFPACK_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${IFPACK_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${IFPACK_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_IFPACK YES)
  find_package_handle_standard_args(IFPACK DEFAULT_MSG IFPACK_LIBRARY)
ENDIF(IFPACK_INCLUDE_PATH AND IFPACK_LIBRARY)

IF(LOCA_INCLUDE_PATH AND LOCA_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${LOCA_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${LOCA_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_LOCA YES)
  find_package_handle_standard_args(LOCA DEFAULT_MSG LOCA_LIBRARY)
ENDIF(LOCA_INCLUDE_PATH AND LOCA_LIBRARY)

IF(ML_INCLUDE_PATH AND ML_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${ML_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${ML_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_ML YES)
  find_package_handle_standard_args(ML DEFAULT_MSG ML_LIBRARY)
ENDIF(ML_INCLUDE_PATH AND ML_LIBRARY)

IF(NOX_INCLUDE_PATH AND NOX_LIBRARY AND NOX_EPETRA_INCLUDE_PATH AND NOX_EPETRA_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${NOX_INCLUDE_PATH} ${NOX_EPETRA_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${NOX_LIBRARY} ${NOX_EPETRA_LIBRARY} ${THYRA_LIBRARY} ${THYRA_EPETRA_LIBRARY}  CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_NOX YES)
  find_package_handle_standard_args(NOX DEFAULT_MSG NOX_LIBRARY)
ENDIF(NOX_INCLUDE_PATH AND NOX_LIBRARY AND NOX_EPETRA_INCLUDE_PATH AND NOX_EPETRA_LIBRARY)

IF(EPETRAEXT_INCLUDE_PATH AND EPETRAEXT_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${EPETRAEXT_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${EPETRAEXT_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_EPETRAEXT YES)
  find_package_handle_standard_args(EPETRAEXT DEFAULT_MSG EPETRAEXT_LIBRARY)
ENDIF(EPETRAEXT_INCLUDE_PATH AND EPETRAEXT_LIBRARY)

IF(KOMPLEX_INCLUDE_PATH AND KOMPLEX_LIBRARY)
  SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${KOMPLEX_INCLUDE_PATH})
  SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${KOMPLEX_LIBRARY} CACHE INTERNAL "Trilinos libraries")
  SET(HAVE_KOMPLEX YES)
  find_package_handle_standard_args(KOMPLEX DEFAULT_MSG KOMPLEX_LIBRARY)
ENDIF(KOMPLEX_INCLUDE_PATH AND KOMPLEX_LIBRARY)

LIST(LENGTH Trilinos_INCLUDE_DIR LEN)
IF(LEN GREATER 1)
    LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIR)
ENDIF(LEN GREATER 1)


#IF(NOT EPETRA_FOUND)
#  MESSAGE(STATUS "Epetra not found.")
#ENDIF()

IF(EPETRA_FOUND AND TEUCHOS_FOUND AND AZTECOO_FOUND)
  SET(Trilinos_FOUND TRUE)
ENDIF()

#IF(KOMPLEX_FOUND)
#ELSE(KOMPLEX_FOUND)
#  MESSAGE(STATUS "Komplex not found.")
#  SET(HAVE_KOMPLEX NO)
#ENDIF()

# Experimental
if(WITH_ZOLTAN)
  IF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_LIBRARY)
    SET(Trilinos_INCLUDE_DIR ${Trilinos_INCLUDE_DIR} ${ZOLTAN_INCLUDE_PATH})
    SET(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} ${ZOLTAN_LIBRARY} CACHE INTERNAL "Trilinos libraries")
    SET(HAVE_ZOLTAN YES)
    find_package_handle_standard_args(ZOLTAN DEFAULT_MSG ZOLTAN_LIBRARY)
  ENDIF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_LIBRARY)
endif(WITH_ZOLTAN)  

IF(Trilinos_FOUND)
  MESSAGE(STATUS "Trilinos packages found.")
#  message("here is LIBS: ${Trilinos_LIBRARIES}")
#  mark_as_advanced (Trilinos_LIBRARIES)
ELSE (Trilinos_FOUND)
  MESSAGE(STATUS "Trilinos packages NOT found in the system (define Trilinos_DIR if installed)")
  IF (Trilinos_FIND_REQUIRED)
    MESSAGE(  FATAL_ERROR 
         "Could not find Trilinos or one of its packages. 
         Please install it according to instructions at\n
         <http://hpfem.org/hermes/doc/src/installation/matrix_solvers/trilinos.html>." )
  ENDIF (Trilinos_FIND_REQUIRED)
ENDIF(Trilinos_FOUND)
