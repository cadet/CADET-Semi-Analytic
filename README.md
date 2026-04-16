[![DOI](https://zenodo.org/badge/38963550.svg)](https://doi.org/10.5281/zenodo.14066869)

# CADET - Semi Analytic Extension

This is the semi-analytic extension of the [CADET project](https://github.com/cadet).
If you find it useful for your own work, we would appreciate acknowledgements of this software and citations of our papers:

* S. Leweke, E. von Lieres (2016): [Fast arbitrary order moments and arbitrary precision solution of the general rate model of column liquid chromatography with linear isotherm](https://doi.org/10.1016/j.compchemeng.2015.09.009), Computers & Chemical Engineering, 84, 350–362.

Please note that the results from the referenced publication can be reproduced exactly using version 1 of the software, which can be [accessed on zenodo](https://doi.org/10.5281/zenodo.14066870).

## Features

* Family of chromatography models: [2DGRM GRM LRMP LRM CSTR](https://cadet.github.io/master/modelling/index.html)) with one component and no or linear (kinetic and rapid-equilibrium) binding
* Arbitrary inlet profiles via piecewise cubic polynomials
* Arbitrary precision solution of the model using a numerical inverse Laplace transform (can optionally be combined with extrapolation)
* Proven error bounds when applicable (see Notes on Error Estimation](doc/interface.md))
* Shared memory parallelization using OpenMP
* Supports the CADET input format, XML and HDF5 as input file formats, CSV for arbitrary precision output
* Command line interface
* Multi-platform: Works on Windows, Linux, and Mac OS X

## Get and build CADET-semi-analytic

CADET-semi-analytic has been successfully built and run on the following platforms:

* Windows
* Linux
* Mac OS X

### Dependencies

* A C++11 capable compiler with OpenMP support (e.g., GCC >= 4.7, Clang >= 3.7, MS Visual Studio >= 2010)
* [CMake](http://cmake.org/)
* [HDF5](http://www.hdfgroup.org/HDF5/)
* [GMP](https://gmplib.org/)
* [MPFR](http://www.mpfr.org/)
* [MPC](http://www.multiprecision.org/)
* [MPFR C++](http://www.holoborodko.com/pavel/mpfr) (header-only, distributed along with CADET-semi-analytic)
* [quadpack++](http://sourceforge.net/projects/quadpackpp) (header-only, distributed along with CADET-semi-analytic)
* [CppAD](http://www.coin-or.org/CppAD/) (header-only, distributed along with CADET-semi-analytic)
* [TCLAP: Templatized C++ Command Line Parser](http://sourceforge.net/projects/tclap/) (header-only, distributed along with CADET-semi-analytic)
* [pugixml](http://code.google.com/p/pugixml/) (distributed along with CADET-semi-analytic)
* Optionally: [FADBAD++](http://www.fadbad.com/) (distributed along with CADET-semi-analytic)

### Building

The project is built with CMake and depends on HDF5, GMP, MPFR, MPC and a C++ toolchain. Two easy approaches:

1) Windows (recommended): use `vcpkg` for dependency management and Visual Studio integration

   - Install `vcpkg` (https://github.com/microsoft/vcpkg) and enable Visual Studio integration:

	 ```powershell
	 .\vcpkg\bootstrap-vcpkg.bat
	 .\vcpkg\vcpkg integrate install
	 .\vcpkg\vcpkg install hdf5[gcc,zlib]:x64-windows gmp mpfr mpc tclap pugixml
	 ```

   - Create an out-of-source build directory and run CMake from the Visual Studio Developer PowerShell or the IDE (vcpkg integration will make packages available automatically):

	 ```powershell
	 mkdir build; cd build
	 cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_FADBAD=OFF
	 cmake --build . --config Release --target INSTALL
	 ```

   - To enable `FADBAD++` or other optional components, set `-DUSE_FADBAD=ON` when calling `cmake` and ensure the library is available via vcpkg or on your system.

2) Linux / macOS: use your package manager or system-provided libraries

   - Install dependencies (example for Debian/Ubuntu):

	 ```bash
	 sudo apt-get update
	 sudo apt-get install build-essential cmake libhdf5-dev libgmp-dev libmpfr-dev libmpc-dev libtclap-dev libpugixml-dev
	 ```

   - Create a build directory and run CMake:

	 ```bash
	 mkdir build && cd build
	 cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_FADBAD=OFF
	 make -j$(nproc)
	 sudo make install
	 ```

## Using CADET-semi-analytic

CADET-semi-analytic uses the [CADET file format](https://cadet.github.io/master/interface/index.html).
However, it is limited to linear single component models. Examples can be found in [CADET-Verification](https://github.com/cadet/CADET-Verification), where CASEMA results are used as reference solutions.

***The interface and program options is documented in the [Interface Documentation](doc/interface.md).***

### Example

Example files are provided under `test/data`.
Suppose your model is contained in the HDF5 file `model.h5`, then you can compute the chromatogram via (windows)
  ```
	casema-cli.exe model.h5 -o chromatogram.csv -e 1e-100 -p 250 -P 20 -t 4
  ```
  The results are written to the file `chromatogram.csv`.

Alternatively, you can choose to write the output to the H5 input file, following the CADET file format.
Note, however, that in this mode, the output precision is constrained by the H5 format, which does not support arbitrary precision, and is instead limited by the numerical limits of the C++ double type.
  ```
	casema-cli.exe model.h5 -e 1e-100 -p 250 -P 20 -t 4
  ```

