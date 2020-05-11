cmake_minimum_required(VERSION 2.6)
project(igampap)

enable_language(CXX)

##################################
# project directories
##################################


#PROJECT_BINARY_DIR
#PROJECT_SOURCE_DIR

##################################
# compiler options
##################################

# C++ compiler

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
#set(CMAKE_CXX_COMPILER "/opt/mylib/ompp-0.8.5/bin/kinst-ompp /usr/bin/g++")
#set(CXX "kinst-ompp g++")
#set(CMAKE_CXX_COMPILER "kinst-ompp g++")

#if(CXX_COMPILER_NAME MATCHES "icc*")
#    message ("Using INTEL C++ compiler")
#    set(COPT "-O3 -Wall -DSELF_CHECK -DUNIX")
#endif()

#if(CXX_COMPILER_NAME MATCHES "gcc*")
#    message ("Using GNU C++ compiler")
#    set(COPT "-std=c++11 -O3 -w -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1  -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math")
#endif()

#message ("${CXX_COMPILER_NAME}")
#message ("${CXX}")

set(COPT "-w -std=c++11 -O3 -g -mavx2 -ftree-vectorize -ftree-vectorizer-verbose=1 -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fopenmp -fpermissive -frounding-math")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPT}")


##################################
# include directories
##################################
include_directories(
src/base
/usr/include/c++/7
/usr/include/eigen3
/usr/include/eigen3/unsupported
/usr/include/vtk-6.3
/usr/include/superlu
/usr/include/openmpi
/usr/include/petsc
)



##################################
# lib directories
##################################
link_directories(
/usr/lib/gcc/x86_64-linux-gnu/7
/usr/lib/x86_64-linux-gnu
/usr/lib
/home/chenna/Documents/myCode/bzfem/lib
)



##################################
# source files
##################################


file(GLOB sources1 src/base/*.h)
file(GLOB sources2 src/base/*.cpp)

##################################
##################################


add_executable(cfdsfemmpi       src/cfdsfemmpi.cpp ${sources1} ${sources2})

target_link_libraries(cfdsfemmpi stdc++ m superlu  blas lapack    metis parmetis  petsc  mpi  mpi_cxx  vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)


install(TARGETS  cfdsfemmpi  RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdsfemmpi/bin)


