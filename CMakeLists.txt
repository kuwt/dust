#Specify minimum version of Cmake
cmake_minimum_required(VERSION 2.6)
enable_language (Fortran)
#set( CMAKE_Fortran_COMPILER mpif90 )
#Project name
project(Dust)

# Check for packages
find_package(OpenMP)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C Fortran)

# Add CGNS installation (for the moment could not make it look for it autonomously)
SET(CGNS_INC "/usr/include" CACHE STRING "Include path for cgns headers") 
SET(CGNS_LIB "/usr/lib" CACHE STRING "Path of the cgns library") 


#Set folders
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})


# Compiler options
#SET(CMAKE_Fortran_COMPILER "/home/matteo/LocalUtils/gcc/install/bin/gfortran")

if (OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(SCOREP_ENABLE "Compile with profiling options" OFF)

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
ENDIF (NOT CMAKE_BUILD_TYPE)



IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")
  set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3")
  #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds -std=f2008ts -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds  -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace")
ELSE()
  message(ERROR "Unknown compiler")
ENDIF()

if (SCOREP_ENABLE)
  SET(CMAKE_Fortran_COMPILER "scorep")
  set(CMAKE_Fortran_FLAGS "gfortran ${CMAKE_Fortran_FLAGS_PROFILE} ${SCOREP_COMPILEFLAGS}")
endif()


# +++++++++++++++++++

## Libraries --------
#find_library(PTHREAD_LIB pthread)
#if(NOT PTHREAD_LIB)
#  message(FATAL_ERROR "pthread library not found")
#endif()
## -lgfortran ???
#
#link_directories(/home/davide/Software/BLAS/lapack-3.7.0/lib)
#link_directories(/home/davide/Software/BLAS/xianyi-OpenBLAS-bcfc298/lib)

# +++++++++++++++++++

#add files
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/geo/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/general/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/aeroel/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/linsys/*.f90
                          ${PROJECT_SOURCE_DIR}/src/io/*.f90
                          ${PROJECT_SOURCE_DIR}/src/options/*.f90)

#set the includes
include_directories(${CGNS_INC} ${HDF5_INCLUDE_DIRS})

#set the link path
link_directories(${CGNS_LIB} ${HDF5_LIBRARY_DIRS})

#add executables to be compiled
add_executable(dust  ${PROJECT_SOURCE_DIR}/src/dust.f90 ${sources})
target_link_libraries(dust ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} cgns)

add_executable(dust_pre  ${PROJECT_SOURCE_DIR}/src/dust_pre.f90 ${sources})
target_link_libraries(dust_pre ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} cgns)
