#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

export SPACK_VERSION=1.1.1
export SPACK_PATH=$ALL_CCFRWORK/spack
export SPACK_USER_PREFIX=$ALL_CCFRWORK/gyselalibxx-spack-install-GENOA
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

module purge

if [ ! -d "$SPACK_PATH" ]; then
    git clone --branch "v${SPACK_VERSION}" --depth 1 https://github.com/spack/spack.git $SPACK_PATH
fi
. $SPACK_PATH/share/spack/setup-env.sh

cp "${TOOLCHAIN_ROOT_DIRECTORY}/packages.yaml" "${SPACK_USER_CONFIG_PATH}"

echo "Preparing the Spack environment..."

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-spack-environment repo update
spack --env gyselalibxx-spack-environment concretize --force
spack --env gyselalibxx-spack-environment install --concurrent-packages 4 --jobs 32
