#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

module purge

export SPACK_USER_PREFIX=$ALL_CCFRWORK/spack-user-install
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache
. $ALL_CCFRWORK/spack/share/spack/setup-env.sh

eval -- "$(
    spack \
        --env gyselalibxx-spack-environment \
        load --sh \
        cmake \
        ddc \
        gcc \
        ginkgo \
        googletest \
        kokkos \
        kokkos-kernels \
        kokkos-tools \
        koliop \
        lapack \
        mpi \
        ninja \
        paraconf \
        pdi \
        pdiplugin-decl-hdf5 \
        pdiplugin-decl-netcdf \
        pdiplugin-mpi \
        pdiplugin-pycall \
        pdiplugin-set-value \
        pdiplugin-trace \
        python \
        py-dask \
        py-h5py \
        py-imageio \
        py-matplotlib \
        py-netcdf4 \
        py-numpy \
        py-pyyaml \
        py-scipy \
        py-sympy \
        py-xarray
)"

module load arch/h100 cuda/12.8.0

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"
