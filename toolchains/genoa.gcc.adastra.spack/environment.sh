#!/bin/bash

ENVIRONMENT_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

module purge

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-GENOA/Configuration.${SPACK_USER_VERSION}"
module use "${SPACK_USER_PREFIX}/modules/tcl/linux-rhel8-zen4"

module load \
    cmake \
    gcc/13.2.1.genoa/zen4/libyaml \
    gcc/13.2.1.genoa/zen4/paraconf \
    gcc/13.2.1.genoa/zen4/pdi \
    gcc/13.2.1.genoa/zen4/pdiplugin-decl-hdf5 \
    gcc/13.2.1.genoa/zen4/pdiplugin-mpi \
    gcc/13.2.1.genoa/zen4/pdiplugin-set-value \
    gcc/13.2.1.genoa/zen4/pdiplugin-trace \
    gcc/13.2.1.genoa/zen4/ginkgo \
    gcc/13.2.1.genoa/zen4/ninja \
    gcc/13.2.1.genoa/zen4/py-matplotlib \
    gcc/13.2.1.genoa/zen4/py-xarray \
    gcc/13.2.1.genoa/zen4/py-h5py

module load cpe/24.07
module load craype-x86-genoa
module load PrgEnv-gnu

module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

module list
