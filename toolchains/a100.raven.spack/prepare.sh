#!/bin/bash


export ALL_CCFRWORK=$HOME
export ALL_CCFRSCRATCH=/ptmp/$USER
mkdir -p $ALL_CCFRSCRATCH/pycache
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache



# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

export SPACK_USER_PREFIX=$ALL_CCFRWORK/spack-user-install
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

module purge

git clone --branch v1.1.0 --depth 1 https://github.com/spack/spack.git $ALL_CCFRWORK/spack || true
. $ALL_CCFRWORK/spack/share/spack/setup-env.sh

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

# Concretize on the compute node
spack --env gyselalibxx-spack-environment install --concurrent-packages 2 --jobs 8
