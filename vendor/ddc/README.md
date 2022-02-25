# The discrete domain computation library (DDC)

## Documentation

See https://maison-de-la-simulation.github.io/ddc/

If you like the project, please leave us a github star.

If you want to know more, join un on [Slack](https://join.slack.com/t/ddc-lib/shared_invite/zt-14b6rjcrn-AwSfM6_arEamAKk_VgQPhg)


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
* PDI
* a c++-17 compiler:
  * gcc-9 should work
