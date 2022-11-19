CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

#################################################
# project information
#################################################
PROJECT(Fortran_COUPLING Fortran CXX C)
STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)
STRING(TOUPPER "Laurent Dupuy" PROJECT_AUTHOR)
                                                                         
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

###################                                                  
# NUMODIS WRAPPER #
###################

# make sure NUMODIS_DIR is defined...
if(NOT DEFINED NUMODIS_DIR)     # ... either as -DNUMODIS_DIR=path_to_numodis
  if(DEFINED ENV{NUMODIS_DIR})  # ... or as an environment variable
    set(NUMODIS_DIR $ENV{NUMODIS_DIR} )
  else()
    message(FATAL_ERROR "=> NUMODIS_DIR must be defined") # ... stop otherwise
  endif()  
endif()

# build numodis and numodis wrapper
add_subdirectory(${NUMODIS_DIR} numodis EXCLUDE_FROM_ALL)
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/numodis )

#################################################
# build driver
#################################################
#file(
#  GLOB
#  driver_files
#  ${PROJECT_SOURCE_DIR}/src/driver.F90
#  )

#add_executable( driver.exe ${driver_files} )

#add_dependencies(driver.exe
#  numodis_fortran_wrapper )

#target_link_libraries( driver.exe
#  numodis_fortran_wrapper
#  )

#set_target_properties( driver.exe
#  PROPERTIES LINKER_LANGUAGE Fortran
#  )

#################################################
# build dynamic library for node export
#################################################
file(
  GLOB
  JavierFiles
  ${PROJECT_SOURCE_DIR}/src/global.F90
  ${PROJECT_SOURCE_DIR}/src/NumodisExportBNodes.F90
  ${PROJECT_SOURCE_DIR}/src/NumodisImport.F90
  ${PROJECT_SOURCE_DIR}/src/NumodisExportStress.F90  
  )

add_library( NumodisCoupling SHARED ${JavierFiles} )

add_dependencies( NumodisCoupling 
  numodis_fortran_wrapper )

SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--disable-new-dtags")

# make sure ELMER_DIR is defined...
if(NOT DEFINED ELMER_DIR)     # ... either as -DELMER_DIR=path_to_elmer
  if(DEFINED ENV{ELMER_HOME})  # ... or as an environment variable
    set(ELMER_DIR $ENV{ELMER_HOME} )
  else()
    message(FATAL_ERROR "=> ELMER_DIR must be defined") # ... stop otherwise
  endif()  
endif()

add_library( elmersolver SHARED IMPORTED )
set_target_properties( elmersolver
  PROPERTIES
  IMPORTED_LOCATION ${ELMER_DIR}/lib/elmersolver/libelmersolver.so )

target_include_directories( NumodisCoupling
  PUBLIC
  ${ELMER_DIR}/share/elmersolver/include )

target_link_libraries( NumodisCoupling
  PUBLIC numodis_fortran_wrapper
  INTERFACE elmersolver
  )

set_target_properties( NumodisCoupling
  PROPERTIES LINKER_LANGUAGE Fortran
  )
