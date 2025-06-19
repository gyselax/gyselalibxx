#!/bin/bash

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-GENOA/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge

module load "${SPACK_USER_VERSION}"
which spack
spack debug report
# Spack must work in a clean, purged environment so it can load modules without
# having to purge itself (which it does not do..). When we spack env activate,
# the same constraint applies.
eval -- "$(spack env activate --prompt --sh gyselalibxx-spack-environment)"

module load cmake
module load cpe/24.07
module load craype-x86-genoa

module load PrgEnv-gnu

module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

module list

unalias despacktivate
unset despacktivate
function despacktivate() {
    eval "$(spack env deactivate --sh)"
}
