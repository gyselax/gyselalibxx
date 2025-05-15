# Pre-made build settings

## General usage

In order to build Gyselalib++ the necessary steps are:

1. Prepare the environment
2. Source the environment
3. Use the environment to build

For example on Adastra the necessary commands are:

```sh
mkdir -p build && cd build
../toolchains/mi250.hipcc.adastra.spack/prepare.sh
source ../toolchains/mi250.hipcc.adastra.spack/environment.sh
cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchains/mi250.hipcc.adastra.spack/toolchain.cmake
```

## Available systems

This folder contains sub-folders with `prepare.sh`, `environment.sh`, and `toolchain.cmake` files. Each sub-folder is designed to work in a different environment. You should choose the folder which best fits your work environment. The available sub-folders are:

- `a100.leonardo.spack` : Designed for use with the A100 hardware available on Cineca's Leonardo supercomputer.
- `cpu.spack.gyselalibxx_env` : Designed for use on a local computer. Installation is handled by spack.
- `docker.gyselalibxx_env` : Designed for use with a docker container. This is used by the CI. It is also possible to use these toolchains if all the dependencies have been installed locally.
    Should you wish to try Gyselalib++ rapidly, the docker container can be found hosted on the GitHub Container Registry : [`ghcr.io/gyselax/gyselalibxx_env:latest`](https://github.com/gyselax/gyselalibxx/pkgs/container/gyselalibxx_env). The provided `environment.sh` file will drop you into a container with access to Gyselalib++'s files.
- `mi250.hipcc.adastra.spack` : Designed for use with the MI250 hardware available on CINES's Adastra supercomputer.
- `v100.persee` : Designed for use with the V100 hardware available on IRFM's Persee cluster.
- `v100.ruche` : Designed for use with the V100 hardware available on Mésocentre's Ruche cluster.
- `xeon.persee` : Designed for use with the CPU hardware available on IRFM's Persee cluster.

## Toolchains

The [toolchains](https://en.wikipedia.org/wiki/Toolchain) are represented using a toolchain file. It summarises CMake build settings for the appropriate machine environment and hardware.

You can use the toolchains like so:

```sh
mkdir -p build && cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchains/<CONFIG>/toolchain.cmake
```

The toolchains files are found in the `toolchains` folder and have the extension `.cmake`. The `CONFIG` folders are named after the environment where they run. E.g. `mi250.hipcc.adastra.spack` contains the files necessary to run on MI250 GPUs using `hipcc` on Adastra. `v100.persee` contains the files necessary to run on V100 GPUs on Persee.

Additionally the folder `docker.gyselalibxx_env` provides toolchains which will work on most platforms. These toolchains are designed to work with the docker environment.

## Environments

Toolchains may require a non default environment. For that purpose environment scripts are provided. Each `CONFIG` folder in `toolchains` should contain (a potentially empty) environment file named `environment.sh`. The environment file is used like so:

```sh
source toolchains/<CONFIG>/environment.sh
```

The *sourcing* of an environment file will generally happen before the build or usage of the built product.

## Preparing environments

The environment may need preparation (for example if the dependencies are not available on the machine). For that purpose the `CONFIG` folders should contain (a potentially empty) preparation file named `prepare.sh`. The preparation file is used like so:

```sh
toolchains/<CONFIG>/prepare.sh
```

## General notes

- The environment and preparation script should be **working directory agnostic**;
- Ideally, if producing files (such as the installation of a dependency), placing them in the appropriate location such as the WORK file system area often found on HPC machines;
- A way to test if a toolchain is valid is to run the `prepare.gyselalibxx.sh` script;
- HPC machines evolve, are notoriously heterogenous, unstable and buggy, you should document issues and workarounds explored as part of the code of the toolchain.
