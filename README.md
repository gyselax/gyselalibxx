# Gyselalib++

Gyselalib++ is a collection of C++ components for writing gyrokinetic semi-lagrangian codes and
similar as well as a collection of such codes.

## Compilation

to compile Gyselalib++:

```
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/gyselalibxx.git gyselalibxx
cd gyselalibxx
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-Wall -Wno-sign-compare" ..
make
```

## Execution

to run the tests:
```
ctest --output-on-failure
```

Then, just have a look at `tests/landau/growthrate_t0.0to45.0.png`:

![tests/landau/fft/growthrate_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/landau/fft/growthrate_t0.0to45.0.png?job=cmake_tests_Release "Landau damping rate")

and `tests/landau/frequency_t0.0to45.0.png`:

![tests/landau/fft/frequency_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/gyselalibxx/-/jobs/artifacts/main/raw/build/tests/landau/fft/frequency_t0.0to45.0.png?job=cmake_tests_Release "Landau damping frequency")

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
