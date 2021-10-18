# The Domain DeComposition library (DDC)

## Compilation

to compile DDC:

```
mkdir build
cd build
cmake ..
make
```

to run the tests:
```
ctest --output-on-failure
```

## Dependencies

* cmake-3.15+
* mdspan
* a c++-17 compiler:
  * gcc-9 should work
