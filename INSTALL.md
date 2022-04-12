
## REQUIREMENTS:

to build dust are required:
- A Fortran compiler
- CMake
- A Lapack/BLAS implementation
- An HDF5 library
- A CGNS library

A Fortran compiler, CMake and a Lapack implementation can be found pre-packed 
in most linux distributions.

## Ubuntu 20.04

#### Compilers
  ```bash
  $ sudo apt install gcc g++ gfortran
  ```

#### Libraries
  ```bash
  $ sudo apt make liblapack-dev libblas-dev libcgns-dev libhdf5-dev libopenblas0
  ```

## DUST building and installation (tested under Ubuntu20.04):

- Create a build folder inside this folder (can be "build" or anything else) and move into it:

  ```bash
  $ mkdir build && cd build
  ```

- Configure cmake with standard options:

  ```bash
  $ cmake -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DWITH_PRECICE=$WITH_PRECICE ../
  ```
  where:
  - **$CMAKE_BUILD_TYPE** can be **Release** or **Debug**
  - **$WITH_PRECICE** can be **ON** or **OFF**

- Build DUST:

  ```bash
  $ make
  ```

- Install DUST (with root privileges if needed):

  ```bash
  $ sudo make install
  ```
  The default install folder should be /usr/local/bin

  Other install folders can be set by setting

  cmake -D CMAKE_INSTALL_PREFIX=/path/to/install/folder ../
