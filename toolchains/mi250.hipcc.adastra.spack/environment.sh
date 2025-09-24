#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge

module load "${SPACK_USER_VERSION}"
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
        ginkgo \
        googletest \
        kokkos \
        kokkos-fft \
        kokkos-kernels \
        kokkos-tools \
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

module load cpe/24.07
module load craype-x86-trento craype-accel-amd-gfx90a
module load PrgEnv-gnu-amd amd-mixed/6.3.3

module list

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"
export MPICH_GPU_SUPPORT_ENABLED=1
