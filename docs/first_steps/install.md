# Gyselalib++

Gyselalib++ is a collection of C++ components for writing gyrokinetic semi-lagrangian codes.
If this is your first interaction with gyselalib++ please check out the documentation on [Getting Started with Gyselalib++](getting_started.md).

## Set-up

In order to set up Gyselalib++ on a new machine, simply run:

```bash
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

or

```bash
git clone --recurse-submodules https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
./bin/install-hooks
```

on a machine for which Gyselalib++ is already used an environment script may be available to set up the necessary modules etc.

Please check the `toolchains/` folder to find the existing configurations. See the documentation about [toolchains](../../toolchains/README.md) for more information on the provided files.

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

The toolchains are found in the folder `toolchains/`. Each toolchain is designed for specific hardware with the exception of those in the folder `docker.gyselalibxx_env/`. These toolchains can be used with the docker environment on most hardware.

For more details about toolchains, see the documentation about [toolchains](../../toolchains/README.md).

## Dependencies

To install dependencies through spack, first follow the the 3 first steps of
<https://github.com/pdidev/spack>

Then execute the following:

```sh
spack env create gyselalibxx spack.yaml
spack env activate gyselalibxx
spack concretize --reuse
spack install
```

For example, you can find a Dockerfile installing these dependencies on ubuntu in
`docker/gyselalibxx_env/Dockerfile`.

## Execution

to run the tests:

```bash
ctest --output-on-failure
```

Then, just have a look at `tests/geometryXVx/landau/fft/growthrate_t0.0to45.0.png`:

![tests/geometryXVx/landau/fft/growthrate\_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/geometryXVx/landau/fft/growthrate_t0.0to45.0.png?job=cmake_tests_Release "Landau damping rate")

and `tests/geometryXVx/landau/fft/frequency_t0.0to45.0.png`:

![tests/geometryXVx/landau/fft/frequency\_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/geometryXVx/landau/fft/frequency_t0.0to45.0.png?job=cmake_tests_Release "Landau damping frequency")
