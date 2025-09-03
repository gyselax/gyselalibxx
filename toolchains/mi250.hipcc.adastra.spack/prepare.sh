#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge
module load "${SPACK_USER_VERSION}"
module list

which spack
spack debug report

mkdir -p -- "${SPACK_USER_CONFIG_PATH}/external-repositories"

# Inject recent PDI recipes into our repository.
git clone https://github.com/pdidev/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" || true
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" fetch
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" checkout 5483cfea7d2d39d654c2962114248f597b3ecf46

git clone https://github.com/gyselax/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" || true
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" fetch
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" checkout releases/v0.23
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" pull

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-spack-environment repo add "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi"
spack --env gyselalibxx-spack-environment repo add "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-spack-environment concretize --force

for ((i = 0; i < 8; ++i)); do
    spack --env gyselalibxx-spack-environment install --jobs 48 &
done

wait

# spack --env gyselalibxx-spack-environment module tcl refresh --delete-tree --yes-to-all
