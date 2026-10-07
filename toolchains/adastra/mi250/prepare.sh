#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

SPACK_USER_VERSION="spack-user-5.0.0"
export SPACK_USER_PREFIX="${ALL_CCFRWORK}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

module purge
module load develop "${SPACK_USER_VERSION}"
module list

which spack
spack debug report

cp "${TOOLCHAIN_ROOT_DIRECTORY}/packages.yaml" "${SPACK_USER_CONFIG_PATH}"

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-spack-environment concretize --force
spack --env gyselalibxx-spack-environment install --concurrent-packages 4 --jobs 32
