#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

WORK_DIRECTORY="/scratch/project_465002992" # The scratch is your typical disk based work.
INSTALL_VERSION="0"
export SPACK_USER_PREFIX="${WORK_DIRECTORY}/gyselalibxx-spack-install-MI250/Configuration.${INSTALL_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/Cache"
# export SPACK_SYSTEM_CONFIG_PATH="${TOOLCHAIN_ROOT_DIRECTORY}/SpackConfiguration"
export SPACK_USER_CONFIG_PATH="${TOOLCHAIN_ROOT_DIRECTORY}/SpackConfiguration"

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX="/flash/project_465002992/pycache"
mkdir -p -- "${PYTHONPYCACHEPREFIX}"

module --force purge
module list

git clone --depth 1 --branch v1.1.1 https://github.com/spack/spack.git "${SPACK_USER_PREFIX}/spack" || true
source "${SPACK_USER_PREFIX}/spack/share/spack/setup-env.sh"

which spack
spack debug report

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-spack-environment concretize --force
spack --env gyselalibxx-spack-environment install --concurrent-packages 4 --jobs 32
