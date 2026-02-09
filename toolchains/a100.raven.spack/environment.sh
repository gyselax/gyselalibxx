#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

export SPACK_USER_PREFIX=$HOME/spack-user-install
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=/ptmp/$USER/pycache

module purge
. $HOME/spack/share/spack/setup-env.sh

eval -- "$(
    spack \
        --env gyselalibxx-spack-environment \
        load --sh \
        cmake \
        ginkgo \
        googletest \
        kokkos \
        kokkos-fft \
        kokkos-kernels \
        kokkos-tools \
        lapack \
        mpi \
        ninja \
        paraconf \
        pdi \
        pdiplugin-decl-hdf5 \
        pdiplugin-decl-netcdf \
        pdiplugin-mpi \
        pdiplugin-set-value \
        pdiplugin-trace \
        python \
        py-dask \
        py-h5py \
        py-imageio \
        py-matplotlib \
        py-netcdf4 \
        py-numpy \
        py-scipy \
        py-sympy \
        py-xarray \
        py-pyyaml
)"

module load gcc/14

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"

export GYSELALIBXX_OPENBLAS_ROOT="$(spack location -i openblas)"
