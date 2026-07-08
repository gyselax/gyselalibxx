#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

WORK_DIRECTORY="/scratch/project_465002992" # The scratch is your typical disk based work.
INSTALL_VERSION="0"
export SPACK_USER_PREFIX="${WORK_DIRECTORY}/gyselalibxx-spack-install-MI250/Configuration.${INSTALL_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/Cache"
# export SPACK_SYSTEM_CONFIG_PATH="${TOOLCHAIN_ROOT_DIRECTORY}/SpackConfiguration"
export SPACK_USER_CONFIG_PATH="${TOOLCHAIN_ROOT_DIRECTORY}/SpackConfiguration"

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX="/flash/project_465002992/pycache"
mkdir -p -- "${PYTHONPYCACHEPREFIX}"

module --force purge
module list

source "${SPACK_USER_PREFIX}/spack/share/spack/setup-env.sh"

which spack
spack debug report

eval -- "$(
    spack \
        --env gyselalibxx-spack-environment \
        load --sh \
        cmake \
        ddc \
        gcc \
        ginkgo \
        gmgpolar \
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

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"
