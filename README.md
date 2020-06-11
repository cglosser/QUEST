[![Build Status](https://travis-ci.org/tbertus/QuEST.svg?branch=development)](https://travis-ci.org/tbertus/QuEST/branch/development) [![DOI](https://zenodo.org/badge/37320094.svg)](https://zenodo.org/badge/latestdoi/37320094)

[![Article](https://img.shields.io/badge/article-Phys.%20Rev.%20A%2096%2C%20033816-B10DC9.svg)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.033816)

# QuEST: The Quantum Electromagnetics Simulation Toolbox
![tube](https://user-images.githubusercontent.com/193990/31498525-20bb0fb0-af30-11e7-9e6a-9a85f99cdb09.png)

Simulation software to model the evolution & electromagnetic interactions of
two-level quantum dots. Designed and built at Michigan State University.

## Prerequisites

* C++14-compatible compiler (tested with [GCC](https://gcc.gnu.org/) and
  [Clang](https://clang.llvm.org/))
* [Eigen3](http://eigen.tuxfamily.org) (at least v3.2.2)
* Boost (at least v.1.55.0):
  * [Program
    Options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)
  * [MultiArray](http://www.boost.org/doc/libs/1_55_0/libs/multi_array/doc/index.html)
  * [Test](http://www.boost.org/doc/libs/1_64_0/libs/test/doc/html/index.html)
* [SILO](https://wci.llnl.gov/simulation/computer-codes/silo) (optional)
* [FFTW](http://www.fftw.org/)

## Building

QuEST relies on [CMake](https://cmake.org/) to generate appropriate compile
scripts. To build the executable, first run

    mkdir build
    cd build

followed by

    cmake $PATH_TO_QUEST && make

This will attempt to build

* `quest` (simulation executable)
* `qtest` (unit test executable)
* `point_gen` (utility to quickly generate distributions of points)
* `siloify` (utility to convert `quest` output to the
  [SILO](https://wci.llnl.gov/simulation/computer-codes/silo) file format for
  use in e.g. [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/))

## Running

QuEST requires three files to specify the parameters of a simulation:
`input.cfg`, `dots.cfg`, and `pulse.cfg`. Of these, `input.cfg` has the most
flexibility; you can specify alternate paths to the other configuration files
as well as modify the simulation parameters (speed of light, number of
particles, timestep, etc.). Run `./quest --help` for details of the possible
parameter options.

`dots.cfg` contains a list of quantum dots, one-per-line, each with the
following format:

```
x y z omega_0 T1 T2 dx dy dz
──┬── ───┬─── ──┬── ───┬────
  │      │      │      └──── transition dipole moment
  │      │      └─────────── decay time constants   
  │      └────────────────── transition frequency
  └───────────────────────── spatial coordinates
```

`pulse.cfg` specifies the incident Gaussian pulse(s) with the following format:

```
E_0 delay sigma omega_L kx ky kz px py pz
─┬─ ──┬── ──┬── ───┬─── ───┬──── ───┬────
 │    │     │      │       │        └──── polarization vector (normalized)
 │    │     │      │       └───────────── wavevector
 │    │     │      └───────────────────── laser frequency
 │    │     └──────────────────────────── pulse width (dimensionless)
 │    └────────────────────────────────── peak shift
 └─────────────────────────────────────── amplitude
```


With all three input files in place, simply run the simulation with `./quest`.
The executable will read everything in, perform the calculation with a
percentage complete indicator, and then produce `output.dat` which contains the
trajectory of the matrix elements for every particle in the system.

## Contributing

Please see CONTRIBUTING.md for details on submitting changes.
