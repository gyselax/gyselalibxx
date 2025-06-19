#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-GENOA/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge
module load "${SPACK_USER_VERSION}"

module load cpe/24.07
module load craype-x86-genoa
module load PrgEnv-gnu

module load cray-fftw
module load cray-python

module list

which spack
spack debug report

eval "$(spack env activate --prompt --sh gyselalibxx-spack-environment)"

unalias despacktivate
unset despacktivate
function despacktivate() {
    eval "$(spack env deactivate --sh)"
}
