## Compilation

to compile voice++:

```
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/voicexx.git
cd voicexx
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Execution

to run the tests:
```
ctest --output-on-failure
```

Then, just have a look at `tests/landau/dampingrate_t0.0to45.0.png`:

![tests/landau/dampingrate_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/voicexx/-/jobs/artifacts/main/raw/build/tests/landau/dampingrate_t0.0to45.0.png?job=cmake_tests_Release "Landau damping rate")

and `tests/landau/frequency_t0.0to45.0.png`:

![tests/landau/frequency_t0.0to45.0.png](https://gitlab.maisondelasimulation.fr/gysela-developpers/voicexx/-/jobs/artifacts/main/raw/build/tests/landau/frequency_t0.0to45.0.png?job=cmake_tests_Release "Landau damping frequency")

## Dependencies

* cmake-3.15+
* a c++-17 compiler:
  * gcc-9+ should work
* Lapack
* PDI
* FFTW-3
* mdspan (embedded)
* googletest (embedded)

For example, you can find a Dockerfile installing these dependencies on ubuntu in 
`voicexx_env/Dockerfile`.

# Architecture

## Operation classes

### `splines` & `boundary_conditions`

Spline interpolation (mesh to spline) & evaluation (spline to mesh)

### `time`

Time iteration, including main time loop & predictor-corrector.

### `vlasov`

Vlasov solver including splitting & advections.

### `efield`

TODO field solver.


