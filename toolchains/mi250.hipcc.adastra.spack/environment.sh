#!/bin/bash

ENVIRONMENT_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

module purge

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
module use "${SPACK_USER_PREFIX}/modules/tcl/linux-rhel8-zen3"

module load \
    cmake \
    gcc/13.2.1.mi250/zen3/libyaml \
    gcc/13.2.1.mi250/zen3/paraconf \
    gcc/13.2.1.mi250/zen3/pdi \
    gcc/13.2.1.mi250/zen3/pdiplugin-decl-hdf5 \
    gcc/13.2.1.mi250/zen3/pdiplugin-mpi \
    gcc/13.2.1.mi250/zen3/pdiplugin-set-value \
    gcc/13.2.1.mi250/zen3/pdiplugin-trace \
    gcc/13.2.1.mi250/zen3/ginkgo \
    gcc/13.2.1.mi250/zen3/eigen \
    gcc/13.2.1.mi250/zen3/ninja \
    gcc/13.2.1.mi250/zen3/py-matplotlib \
    gcc/13.2.1.mi250/zen3/py-xarray \
    gcc/13.2.1.mi250/zen3/py-h5py

unset HIPCC_COMPILE_FLAGS_APPEND

module load cpe/24.07
# FIXME:
# craype-accel-amd-gfx90a Error with the cray wrappers clang: error: unsupported option '-fopenmp-targets=' for language mode 'HIP'
module load craype-x86-trento
# NOTE: Force 6.3.3 due to startup failures (https://github.com/gyselax/gyselalibxx/pull/198#issuecomment-2943081411)
module load PrgEnv-gnu-amd amd-mixed/6.3.3
module load rocm/6.3.3
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

module list
