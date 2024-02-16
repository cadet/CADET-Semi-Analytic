![CADET Logo](doc/logo/CADET-GitHub.png "Chromatography Analysis and Design Toolkit")

# Chromatography Analysis and Design Toolkit - Semi Analytic Extension

The Chromatography Analysis and Design Toolkit (CADET) is developed at the Institute of Bio- and Geosciences 1 (IBG-1) of Forschungszentrum Jülich (FZJ) under supervision of Dr. Eric von Lieres. This is the semi analytic extension of the [core CADET project](https://github.com/modsim/cadet), of which both are freely distributed (under the terms of the GPLv3) as a contribution to the scientific community. If you find it useful for your own work, we would appreciate acknowledgements of this software and citations of our papers:

* S. Leweke, E. von Lieres (2016): [Fast arbitrary order moments and arbitrary precision solution of the general rate model of column liquid chromatography with linear isotherm](https://doi.org/10.1016/j.compchemeng.2015.09.009), Computers & Chemical Engineering, 84, 350–362.

## Features

* Fast arbitrary order moments using the Laplace transform of the GRM, algorithmic differentiation, and extrapolation
* Arbitrary inlet profiles via piecewise cubic polynomials in moment calculation
* Arbitrary precision solution of the GRM using a numerical inverse Laplace transform (can optionally be combined with extrapolation)
* Proven error bounds for GRM solutions in case of quasi-stationary binding
* Processes one component general rate models (GRM) with linear isotherm
* Suited for dynamic and quasi-stationary binding
* Shared memory parallelization using OpenMP
* Supports XML and HDF5 as input formats, CSV for output
* Command line interface
* Multi-platform: Works on Windows, Linux, and Mac OS X

## Get and build CADET-semi-analytic

CADET-semi-analytic has been successfully built and run on the following platforms:

* Linux (Ubuntu 14.04, SLES11 SP3)
* Mac OS X (10.6.8 Snow Leopard, 10.9 Mavericks).

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

1. Download and build the requirements. Note that CADET-semi-analytic needs the HDF5-C++ library which is not built by default (use `--enable-cxx` when calling `./configure`).
2. Get the [latest source code](https://github.com/modsim/cadet-semi-analytic/archive/master.zip) of CADET-semi-analytic
3. Unpack to the folder `casema/code`
4. Create a new folder `casema/build` and change to it
5. Call `cmake` and use the environment variables `GMP_ROOT`, `MPFR_ROOT`, `MPC_ROOT`, and `HDF5_ROOT` to point CMake to the top level directories of the installed libraries if it does not find them automatically. Use `-DCMAKE_INSTALL_PREFIX` to tell CMake the installation directory and `DCMAKE_CXX_COMPILER` if you want to use a non-default compiler. On Linux a suitable command might look like this (paths need to be adjusted):

  ```
GMP_ROOT=libs/gmp/ MPFR_ROOT=libs/mpfr/ MPC_ROOT=libs/mpc/ HDF5_ROOT=libs/hdf5/ cmake -DCMAKE_INSTALL_PREFIX=install/ -DCMAKE_CXX_COMPILER=g++-4.8 -DCMAKE_C_COMPILER=gcc-4.8 ../code
  ```
  Note that you need to pass the switch `-DUSE_FADBAD=ON` in order to enable FADBAD++ support.
6. After CMake has successfully run, the software is compiled and installed by executing `make install`

## Using CADET-semi-analytic

CADET-semi-analytic uses the same file format [CADET](https://github.com/modsim/cadet) does. However, it is limited to one component and linear isotherms. Example models taken from S. Qamar et al. (2014) "Analytical solutions and moment analysis of general rate model for linear liquid chromatography" (Chemical Engineering Science, 107, 192–205, doi: [10.1016/j.ces.2013.12.019](http://dx.doi.org/10.1016/j.ces.2013.12.019)) are provided in the folder `examples`.

Suppose your model is contained in the HDF5 file `model.h5`, then you can do the following things:

* Compute the chromatogram via

  ```
chrom -i model.h5 -o chromatogram.csv -e 1e-100 -p 250 -P 20 -t 4 --convabs 0 --convrel 0
  ```
  where no extrapolation method is used and, hence, the convergence detection tolerances are set to 0. This command also requests the usage of 250 decimal digits precision arithmetics (but only 20 digits of them are written to file), parallelization using 4 threads, and the total error to be less than 10^(-100). Extrapolation is enabled by adding `-x MET` to the command line, where `MET` is one of `ide`, `ads`, `wem`, `wrm`, `iad`, `lum`, `ltm`, `ibt`, `btm`, `nam`, `rem`, or `sgr`. The results are written to the file `chromatogram.csv`.

* Compute moments via Laplace transform by calling

  ```
moments -i model.h5 -o moments.csv -p 5000 -P 20 --convabs 1e-100 --convrel 1e-100 -m 100 -n 100 -v --stopsingleconv -e wem
  ```
  where the precision is set to 5000 decimal digits (of which 20 are output), the convergence tolerances are set to 10^(-100), Wynn's epsilon algorithm is used for extrapolation, and a maximum of 100 iterations are performed. The extrapolation stops after convergence has been detected, thanks to the `--stopsingleconv` flag and all intermediate results are written to file because of the `-v` flag. This command computes the first 100 moments (`-m 100`) and writes them to the file `moments.csv`.

* Calculate analytical moments up to fourth order for quasi-stationary binding from formulas derived by Qamar et al. (2014):

  ```
moments -i model.h5 -o moments.csv -p 5000 -P 20 -a
  ```
  where the flags have the same meaning as before. **Attention:** A comparison indicates that the analytical formulas for the third and fourth non-central as well as the fourth central moment are impaired, presumably by typographical errors!

* Calculate moments via numerical inverse Laplace transform and adaptive Gauss-Kronrod quadrature using

  ```
momentQuadrature -i model.h5 -o moments.csv -p 5000 -P 20 -e 1e-100 --quadrel 1e-50 --quadabs 1e-50 -N 10 -M 50 -m 30
  ```
  where the precision flags `-P`, `-p` are the same as before, the inverse Laplace transform is calculated with at most 1^(-100) error, the adaptive integration tries to meet 1^(-50) absolute and relative error tolerances (`--quadrel`, `--quadabs`) and may use up to 50 (`-M`) subdivisions and a Gauss-Kronrod rule comprising 61 = 2*30+1 points (`-m`). The command computes the 10th (`-N`) central and non-central moments and writes them to the file `moments.csv`.
