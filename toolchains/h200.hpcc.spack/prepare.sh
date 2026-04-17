#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

export SPACK_USER_PREFIX=$HOME/spack-user-install
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=/ptmp/$USER/pycache

module purge
module load gnu/gcc-12.3
module load openmpi/gpu/4.1.5
module load cuda/12.9

git clone --branch v1.1.0 --depth 1 https://github.com/spack/spack.git $HOME/spack || true
. $HOME/spack/share/spack/setup-env.sh

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

# Concretize on the compute node
spack --env gyselalibxx-spack-environment install --concurrent-packages 4 --jobs 16
