# Gyselalib++

Gyselalib++ is a collection of C++ components for writing gyrokinetic semi-lagrangian codes and
similar as well as a collection of such codes.
It is based on [DDC](https://ddc.mdls.fr/). We strongly encourage new developers to begin by reading our documentation about [Using DDC in Gyselalibxx](./docs/DDC_in_gyselalibxx.md).

## Compilation

to compile Gyselalib++:

```
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=<TOOLCHAIN_FILE> ..
make
```

The toolchains are found in the folder `toolchains/`. Each toolchain is designed for specific hardware with the exception of those in the folder `docker.gyselalibxx_env/`. These toolchains can be used with the docker environment on most hardware.

For more details about toolchains, see the documentation about [toolchains](./toolchains/README.md).

## Execution

to run the tests:
```
ctest --output-on-failure
```

Then, just have a look at `tests/geometryXVx/landau/fft/growthrate_t0.0to45.0.png`:

![tests/geometryXVx/landau/fft/growthrate\_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/geometryXVx/landau/fft/growthrate_t0.0to45.0.png?job=cmake_tests_Release "Landau damping rate")

and `tests/geometryXVx/landau/fft/frequency_t0.0to45.0.png`:

![tests/geometryXVx/landau/fft/frequency\_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/geometryXVx/landau/fft/frequency_t0.0to45.0.png?job=cmake_tests_Release "Landau damping frequency")

## Dependencies

To install dependencies through spack, first follow the the 3 first steps of 
https://github.com/pdidev/spack

Then execute the following:
```sh
spack env create gyselalibxx spack.yaml
spack env activate gyselalibxx
spack concretize --reuse
spack install
```

For example, you can find a Dockerfile installing these dependencies on ubuntu in 
`docker/gyselalibxx_env/Dockerfile`.
