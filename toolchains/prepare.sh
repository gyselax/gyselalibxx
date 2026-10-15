#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    exit 1
fi

if [[ $# -lt 3 ]]; then
    echo "Usage: bash prepare.sh <project-name> <platform-name> <arch-name>" >&2
    exit 1
fi
GYSELA_PROJECT_NAME="$1"
GYSELA_PLATFORM_NAME="$2"
GYSELA_ARCH_NAME="$3"

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

. ${TOOLCHAIN_ROOT_DIRECTORY}/platforms/${GYSELA_PLATFORM_NAME}/config.sh

GYSELA_SPACK_PREFIX=${GYSELA_SPACK_BASE}/gysela-spack-${GYSELA_SPACK_VERSION}
if [ ! -d "${GYSELA_SPACK_PREFIX}" ]; then
    mkdir --parents "${GYSELA_SPACK_PREFIX}"
    chmod g+s "${GYSELA_SPACK_PREFIX}"
    git clone --branch v${GYSELA_SPACK_VERSION} --depth 1 https://github.com/spack/spack.git "${GYSELA_SPACK_PREFIX}"
fi

export SPACK_DISABLE_LOCAL_CONFIG=True
export PYTHONDONTWRITEBYTECODE=True
. ${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh

env="${GYSELA_PROJECT_NAME}-${GYSELA_ARCH_NAME}"
echo "Preparing the Spack environment ${env}"

spack env remove --yes-to-all "${env}"
spack env create "${env}"

env_dir="$(spack location --env ${env})"
cp "${TOOLCHAIN_ROOT_DIRECTORY}/projects/${GYSELA_PROJECT_NAME}.yaml" "${env_dir}/spack.yaml"
mkdir --parents "${env_dir}/gysela"
cp "${TOOLCHAIN_ROOT_DIRECTORY}/platforms/${GYSELA_PLATFORM_NAME}/${GYSELA_ARCH_NAME}/externals.yaml" "${env_dir}/"
cp "${TOOLCHAIN_ROOT_DIRECTORY}/platforms/${GYSELA_PLATFORM_NAME}/${GYSELA_ARCH_NAME}/packages.yaml" "${env_dir}/gysela"
cp "${TOOLCHAIN_ROOT_DIRECTORY}/platforms/${GYSELA_PLATFORM_NAME}/${GYSELA_ARCH_NAME}/repos.yaml" "${env_dir}/gysela"

spack env activate "${env}"
spack repo update
spack concretize --quiet
spack spec --install-status --namespaces
spack install --jobs "${GYSELA_SPACK_JOBS}"
spack env deactivate

find "${GYSELA_SPACK_PREFIX}" \( -type d -type f \) -user "${USER}" \
    -exec chgrp "${GYSELA_SPACK_GROUP}" {} + \
    -exec chmod g+rwX {} +
