#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

SPACK_USER_VERSION="spack-user-5.0.0"

export SPACK_USER_PREFIX="${ALL_CCFRWORK}/gyselalibxx-spack-install-GENOA/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

module purge

module load develop "${SPACK_USER_VERSION}"
which spack
spack debug report
# Spack must work in a clean, purged environment so it can load modules without
# having to purge itself or clearing environment variables (which it does not
# do..). When we spack env activate, the same constraint applies.
# Use spack load instead of an environment activation as it should limit the
# inode produced by the environment's view.
# eval -- "$(spack env activate --prompt --sh gyselalibxx-spack-environment)"
# unalias despacktivate
# unset despacktivate
# function despacktivate() {
#     eval "$(spack env deactivate --sh)"
# }

eval -- "$(
    spack \
        --env gyselalibxx-spack-environment \
        load --sh \
        cmake \
        gcc \
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

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"

export GYSELALIBXX_OPENBLAS_ROOT="$(spack location -i openblas)"
