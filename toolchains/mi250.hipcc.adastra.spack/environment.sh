#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
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
module load craype-x86-trento craype-accel-amd-gfx90a
# NOTE: Force 6.3.3 due to startup failures (https://github.com/gyselax/gyselalibxx/pull/198#issuecomment-2943081411)
module load PrgEnv-gnu-amd amd-mixed/6.3.3
module load rocm/6.3.3
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

module list

unalias despacktivate
unset despacktivate
function despacktivate() {
    eval "$(spack env deactivate --sh)"
}
