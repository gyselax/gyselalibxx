# Gyselalib++ Installation

This document provides comprehensive instructions for installing, compiling, and running Gyselalib++ on various systems.

## Table of Contents

1. [Cloning the Repository](#cloning-the-repository)
2. [Environment Setup](#environment-setup)
3. [Compilation](#compilation)
4. [Dependencies Management](#dependencies-management)
5. [Running Tests](#running-tests)

## Cloning the Repository

To get started with Gyselalib++, you'll need to clone the repository along with its submodules. You can do this using either SSH or HTTPS:

### SSH Method (recommended for developers with SSH keys configured)

```bash
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

### HTTPS Method

```bash
git clone --recurse-submodules https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

The `install-hooks` script sets up necessary git hooks for the repository.

## Environment Setup

Gyselalib++ requires specific environment configurations depending on your system. Pre-configured environment scripts are available in the `toolchains/` directory.

### Available Toolchains

The toolchains are organised by hardware architecture. For detailed information about available toolchains, see the [toolchains documentation](../../toolchains/README.md).

On a machine for which Gyselalib++ is already used an environment script may be available to set up the necessary modules etc.

### Example: Adastra Supercomputer Setup

For example in order to set up the environment on the Adastra supercomputer simply run:

```bash
source toolchains/mi250.hipcc.adastra.spack/prepare.sh
source toolchains/mi250.hipcc.adastra.spack/environment.sh
```

## Compilation

To compile Gyselalib++, navigate to the folder where the library was cloned to and run:

```bash
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=<TOOLCHAIN_FILE> ..
make
```

### Docker Toolchains

The toolchains are found in the folder `toolchains/`. Each toolchain is designed for specific hardware with the exception of those in the folder `docker.gyselalibxx_env/`. These toolchains can be used with the docker environment on most hardware.

For more details about toolchains, see the documentation about [toolchains](../../toolchains/README.md).

## Dependencies Management

1. To install dependencies through spack, first follow the the 3 first steps of
<https://github.com/pdidev/spack>
2. Then execute the following:

```sh
spack env create gyselalibxx spack.yaml
spack env activate gyselalibxx
spack concretize --reuse
spack install
```

### Docker Example

For example, you can find a Dockerfile installing these dependencies on Ubuntu at
`docker/gyselalibxx_env/Dockerfile`.

## Running Tests

After successful compilation, you can run the test suite with:

```bash
ctest --output-on-failure
```

### Expected Test Results

Upon successful execution, you should see output plots demonstrating the Landau damping characteristics:

1. **Landau Damping Rate**
    Growth Rate Plot: `tests/geometryXVx/landau/fft/growthrate\_t0.0to45.0.png`
2. **Landau Damping Frequency**
    Frequency Plot: `tests/geometryXVx/landau/fft/frequency\_t0.0to45.0.png`