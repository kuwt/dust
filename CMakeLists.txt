#Specify minimum version of Cmake
cmake_minimum_required(VERSION 2.6)
enable_language (Fortran)
#Project name
project(Dust)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#Set folders
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)


# =========== COMPILER ========================
#SET(CMAKE_Fortran_COMPILER "/home/matteo/LocalUtils/gcc/install/bin/gfortran")

find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -g -fopenmp")
  set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3 -g -fopenmp")
  #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds -std=f2008ts -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ggdb -g -O0 -fopenmp -pedantic -Wall -Wconversion -Wconversion-extra -fcheck=all -fcheck=bounds  -fall-intrinsics -finit-real=snan -ffpe-trap=invalid -finit-integer=-666 -Wintrinsics-std -fbacktrace")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -standard-semantics")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -g -xHost")
  set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O3 -g -xHost")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g -O0 -stand f08 -check all -check noarg_temp_created -warn all -fpe0 -traceback -ftrapuv -fp-stack-check")
ELSE()
  message(ERROR "Unknown compiler")
ENDIF()


#LIBRARIES
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

#Lapack
OPTION(DUST_MKL "Use MKL libraries instead of BLAS/LAPACK" OFF)
IF (DUST_MKL)
#  SET(BLA_VENDOR "Intel") #search only for intel libraries
  find_package(MKL REQUIRED)
  LIST(APPEND incldirs ${MKL_INCLUDE_DIRS}) 
  LIST(APPEND linkedlibs ${MKL_LIBRARIES} pthread dl)
  IF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") 
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGSi_RELEASE} -mkl")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -mkl")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -mkl")
  ENDIF() 
ELSE()
  find_package(LAPACK REQUIRED)
  LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
ENDIF()

#HDF5
find_package(HDF5 REQUIRED COMPONENTS C Fortran)
LIST(APPEND linkedlibs ${HDF5_LIBRARIES})

# Add CGNS installation (for the moment could not make it look for it autonomously)
SET(CGNS_INC "/usr/include" CACHE STRING "Include path for cgns headers") 
SET(CGNS_LIB "/usr/lib" CACHE STRING "Path of the cgns library") 
LIST(APPEND linkedlibs cgns)

#set the includes
include_directories(${incldirs} ${CGNS_INC} ${HDF5_INCLUDE_DIRS})

#set the link path
link_directories(${CGNS_LIB} ${HDF5_LIBRARY_DIRS})



#add files
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/geo/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/general/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/aeroel/*.f90 
                          ${PROJECT_SOURCE_DIR}/src/linsys/*.f90
                          ${PROJECT_SOURCE_DIR}/src/io/*.f90
                          ${PROJECT_SOURCE_DIR}/src/options/*.f90)


#create a library
add_library(dustlibF90 OBJECT ${sources})
add_library(dustlib STATIC $<TARGET_OBJECTS:dustlibF90>)
target_link_libraries(dustlib)

#add executables to be compiled
add_executable(dust  ${PROJECT_SOURCE_DIR}/src/dust.f90)
add_dependencies(dust dustlib)
target_link_libraries(dust dustlib ${linkedlibs})

add_executable(dust_pre  ${PROJECT_SOURCE_DIR}/src/dust_pre.f90)
add_dependencies(dust_pre dustlib)
target_link_libraries(dust_pre  dustlib ${linkedlibs})

add_executable(dust_post  ${PROJECT_SOURCE_DIR}/src/dust_post.f90)
add_dependencies(dust_post dustlib)
target_link_libraries(dust_post  dustlib ${linkedlibs})
