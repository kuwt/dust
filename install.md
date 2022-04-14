
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

<details>
  <summary markdown="span">Compilers</summary>

#### Compilers
  ```bash
  $ sudo apt install gcc g++ gfortran
  ```
</details>

<details>
  <summary markdown="span">Libraries</summary>

#### Libraries
  ```bash
  $ sudo apt install liblapack-dev libblas-dev libcgns-dev libhdf5-dev libopenblas0
  ```
</details>

<details>
  <summary markdown="span">Installation</summary>
  
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
  - **$WITH_PRECICE** can be **YES** or **NO**

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
</details>

## Coupling with preCICE-MBDyn

Compile DUST with **$WITH_PRECICE**=**ON**.

<details>
  <summary markdown="span">preCICE</summary>

#### preCICE
Visit <https://precice.org/quickstart.html>

</details>

<details>
  <summary markdown="span">MBDyn</summary>

#### MBDyn
Visit <https://www.mbdyn.org/?Software_Installation>. 

MBDyn must be compiled on branch **develop** with 
the following configure command:
 ```bash
  $ ./configure --enable-netcdf --with-lapack --enable-python
  ```
</details>

