# Pre-made build settings

## General usage

In order to build Gyselalib++ the necessary steps are:
1. Prepare the environment
2. Source the environment
3. Use the environment to build

For example on Adastra the necessary commands are:

```
$ mkdir -p build && cd build
$ ../toolchains/mi250.hipcc.adastra.spack/prepare.sh
$ source ../toolchains/mi250.hipcc.adastra.spack/environment.sh
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchains/mi250.hipcc.adastra.spack/toolchain.cmake
```

## Toolchains

The [toolchains](https://en.wikipedia.org/wiki/Toolchain) are represented using a toolchain file. It summarises CMake build settings for the appropriate machine environment and hardware.

You can use the toolchains like so:

```
$ mkdir -p build && cd build
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchains/<CONFIG>/toolchain.cmake
```

The toolchains files are found in the `toolchains` folder and have the extension `.cmake`. The `CONFIG` folders are named after the environment where they run. E.g. `mi250.hipcc.adastra.spack` contains the files necessary to run on MI250 GPUs using `hipcc` on Adastra. `v100.persee` contains the files necessary to run on V100 GPUs on Persee.

Additionally the folder `docker.gyselalibxx_env` provides toolchains which will work on most platforms. These toolchains are designed to work with the docker environment.

## Environments

Toolchains may require a non default environment. For that purpose environment scripts are provided. Each `CONFIG` folder in `toolchains` should contain (a potentially empty) environment file named `environment.sh`. The environment file is used like so:

```
$ source toolchains/<CONFIG>/environment.sh
```

The *sourcing* of an environment file will generally happen before the build or usage of the built product.


## Preparing environments

The environment may need preparation (for example if the dependencies are not available on the machine). For that purpose the `CONFIG` folders should contain (a potentially empty) preparation file named `prepare.sh`. The preparation file is used like so:

```
$ toolchains/<CONFIG>/prepare.sh
```

## General notes

- The environment and preparation script should be **working directory agnostic**;
- Ideally, if producing files (such as the installation of a dependency), placing them in the appropriate location such as the WORK file system area often found on HPC machines;
- A way to test if a toolchain is valid is to run the `prepare.gyselalibxx.sh` script;
- HPC machines evolve, are notoriously heterogenous, unstable and buggy, you should document issues and workarounds explored as part of the code of the toolchain.
