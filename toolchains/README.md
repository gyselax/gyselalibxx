# Pre-made build settings

## Toolchains

The [toolchains](https://en.wikipedia.org/wiki/Toolchain) are represented using a toolchain file. It summarize CMake build settings for the appropriate machine environment and hardware.

On can use the toolchains like so:
```
$ cd koliop && mkdir -p build && cd build
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchains/mi250.amd.adastra/toolchain.cmake
```
