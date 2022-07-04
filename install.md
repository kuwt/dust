# Binary packages
Pre-compiled binary packages are available for Ubuntu 20.04 (focal) and 22.04 (jammy). They include only the standalone DUST version, without the preCICE-MBDyn coupling.

For Ubuntu, you can download and install DUST as follows:
- Download the related package from the [release page](https://gitlab.com/dust_project/dust_dev/-/releases/)
- Move to the directory where you download the package and run
  ```bash
    sudo apt install ./dust_0.7.2b-1_amd64_focal.deb
  ```
Change `focal` to `jammy` in the line above according to your version.
# Build DUST from source
If you need the coupled preCICE-MBDyn version of DUST, or if you wish to develop the software, you have to build it from source.
## Requirements:

To build DUST are required:

- a Fortran compiler
- CMake 
- a Lapack/BLAS implementation
- an HDF5 library
<br/>supported and tested versions: <ins>1.10</ins> , <ins> 1.12</ins><br/>
- a CGNS library
<br/>supported and tested versions: <ins>3.4.0</ins> , <ins>4.3.0</ins><br/>



A Fortran compiler, CMake and a Lapack implementation can be found pre-packed 
in most Linux distributions.

## Ubuntu 20.04 or above

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
  $ sudo apt install liblapack-dev libblas-dev libopenblas-dev libopenblas0 libcgns-dev libhdf5-dev
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

- If you want to compile with intel MKL libraries first:
  ```bash
    source /opt/intel/oneapi/setvars.sh
  ```
  then: 
  ```bash
  $ cmake -DDUST_MKL=YES ../ 
  ``` 
  
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

Compile DUST with **$WITH_PRECICE**=**ON** and include the adapter and interface to your Python path.

For example, add these line to your ~/.bashrc file:

  ```bash
  $ export PYTHONPATH="/path/to/dust/utils/adapter":$PYTHONPATH
  ```

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

