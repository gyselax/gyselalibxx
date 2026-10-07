#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

BASE_SPACK=${ALL_CCFRSCRATCH}/gysela-spack
SPACK_VERSION="1.2.2"
export SPACK_PREFIX=${BASE_SPACK}/spack-${SPACK_VERSION}
export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CACHE_PATH=${BASE_SPACK}/cache
export PYTHONDONTWRITEBYTECODE=True

git clone --branch v1.2.2 --depth 1 https://github.com/spack/spack.git "${SPACK_PREFIX}" || true

. ${SPACK_PREFIX}/share/spack/setup-env.sh

echo "Preparing the Spack environments..."

spack env remove --yes-to-all gyselalibxx-mi250
spack env create gyselalibxx-mi250 "${TOOLCHAIN_ROOT_DIRECTORY}/mi250/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-mi250 repo update
spack --env gyselalibxx-mi250 concretize
spack --env gyselalibxx-mi250 install --jobs 96

spack env remove --yes-to-all gyselalibxx-mi300
spack env create gyselalibxx-mi300 "${TOOLCHAIN_ROOT_DIRECTORY}/mi300/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-mi300 repo update
spack --env gyselalibxx-mi300 concretize
spack --env gyselalibxx-mi300 install --jobs 96

spack env remove --yes-to-all gyselalibxx-genoa
spack env create gyselalibxx-genoa "${TOOLCHAIN_ROOT_DIRECTORY}/genoa/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-genoa repo update
spack --env gyselalibxx-genoa concretize
spack --env gyselalibxx-genoa install --jobs 96
