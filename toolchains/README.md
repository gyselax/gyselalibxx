# Gyselalib++ Installation

This document provides comprehensive instructions for installing, compiling, and running Gyselalib++ on various systems.

## Table of Contents

1. [Cloning the Repository](#cloning-the-repository)
2. [Dependencies Management](#dependencies-management)
3. [Environment Setup](#environment-setup)
4. [Compilation](#compilation)
5. [Running Tests](#running-tests)

## Cloning the Repository

To get started with Gyselalib++, you'll need to clone the repository along with its submodules. You can do this using either SSH or HTTPS:

### SSH Method (required for developers who want to push to the repository)

```bash
git clone --recurse-submodules git@github.com:gyselax/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

### HTTPS Method

```bash
git clone --recurse-submodules https://github.com/gyselax/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

The `install-hooks` script sets up necessary git hooks for the repository.

## Dependencies Management

Dependencies can be installed using a script called `prepare.sh`. The best way to install dependencies is machine dependent, so the correct file will depend on the machine where you want to run.

Ideally this script should be run exactly once per machine. The script should be **working directory agnostic**, in other words it should behave in the same way regardless of the folder it is run from. The preparation file is used like so:

```sh
./toolchains/<CONFIG>/prepare.sh
```

:warning: Installation can take a few hours.

:warning: Where possible, on commonly-used machines we install the dependencies in a location that is accessible for all users. In this case no `prepare.sh` script is provided. Only the administrator who controls the location where the dependencies were installed will need to look at installation scripts.

:warning: Ideally, if files are produced by the prepare script (such as the installation of a dependency), these should be placed in the appropriate location for the machine. For example, in the WORK file system area often found on HPC machines.

### Available systems

The folder [toolchains/](https://github.com/gyselax/gyselalibxx/tree/devel/toolchains) contains sub-folders with files that can be used for dependency management (`prepare.sh`, environment setup (`environment.sh`), and compilation `toolchain.cmake`. Each sub-folder is designed to work in a different environment. You should choose the folder which best fits your work environment. The available sub-folders are:

- `a100.leonardo.spack` : Designed for use with the A100 hardware available on Cineca's Leonardo supercomputer.
- `cpu.spack.gyselalibxx_env` : Designed for use on a local computer. Installation is handled by spack.
- `docker.gyselalibxx_env` : Designed for use with a docker container. This is used by the CI.
    Should you wish to try Gyselalib++ rapidly, the docker container can be found hosted on the GitHub Container Registry : [`ghcr.io/gyselax/gyselalibxx_env:latest`](https://github.com/gyselax/gyselalibxx/pkgs/container/gyselalibxx_env). The provided `environment.sh` file will drop you into a container with access to Gyselalib++'s files.
- `h100.jean-zay.spack` : Designed for use with the H100 hardware available on IDRIS' Jean-Zay supercomputer. It installs a new Spack instance in the shared work directory of the project. It may need to be adapted depending on your needs.
- `mi250.hipcc.adastra.spack` : Designed for use with the MI250 hardware available on CINES's Adastra supercomputer.
- `persee/v100` : Designed for use with the V100 hardware available on IRFM's Persee cluster.
- `persee/xeon` : Designed for use with the CPU hardware available on IRFM's Persee cluster.
- `v100.ruche` : Designed for use with the V100 hardware available on Mésocentre's Ruche cluster.

## Environment Setup

Once all the dependencies have been installed the environment must be set up in order to access the installed components. This can be done using a script called `environment.sh`. The steps required to set up depend on how the dependencies were installed, so the correct file will depend on the machine where you want to run.

For this purpose, environment scripts are provided. Each `CONFIG` folder in `toolchains` should contain (a potentially empty) environment file named `environment.sh`. The script should be **working directory agnostic**, in other words it should behave in the same way regardless of the folder it is run from. The environment file is used like so:

```sh
source toolchains/<CONFIG>/environment.sh
```

The *sourcing* of an environment file will generally happen before the build or usage of the built product.

:warning: This script modifies environment variables so it must be executed every time you connect to a cluster, or everytime you open a new terminal on a local machine.

## Compilation

Once the dependencies are installed and the environment has been set up the compilation is mostly machine independent. There are however a few parameters which may depend on the machine. We use [toolchains](https://en.wikipedia.org/wiki/Toolchain) to describe the CMake build settings for the appropriate machine environment and hardware.

You can use the toolchains by running the following command from the project root directory:

```sh
cmake -B build -DCMAKE_TOOLCHAIN_FILE=toolchains/<FOLDER>/<CONFIG>.cmake .
```

Kokkos options are one notable reason for the machine dependency. The correct backend must be activated to benefit from the acceleration available on the machine. In order to ensure that the best optimisation flags are used, a `toolchain.cmake` file is provided for each cluster. When using a local installation, or when a non-optimised compilation is wanted (e.g. a debug compilation) a less specific toolchain can be used. You can find common toolchain choices in the folder `toolchains/common_toolchains`.

:warning: Developers are advised to work with a debug configuration. It is usually simpler to begin with a serial debug configuration and only move onto more complex configurations (e.g. machine specific GPU debug configuration) once the CPU version is working. Serial configurations are found in the folder `toolchains/common_toolchains`.

## Running Tests

After successful compilation, you can run the test suite with:

```bash
ctest --output-on-failure
```

If you want to run a specific test you can use the ctest flag `-R` which will use a regex to search for tests with a name containing the specified keyword. For example, for tests related to splines:

```bash
ctest -R Splines --verbose --output-on-failure
```

### Expected Test Results

Upon successful execution, you should see all tests passing. You can also find output plots demonstrating the Landau damping characteristics:

1. **Landau Damping Rate**
    Growth Rate Plot: `tests/geometryXVx/landau/fft/growthrate\_t0.0to45.0.png`
2. **Landau Damping Frequency**
    Frequency Plot: `tests/geometryXVx/landau/fft/frequency\_t0.0to45.0.png`

## General notes

- HPC machines evolve, are notoriously heterogenous, unstable and buggy, you should document issues and workarounds explored as part of the code of the toolchain.
