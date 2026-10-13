#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

GYSELA_SPACK_GROUP="gen2224"
GYSELA_SPACK_JOBS=96
GYSELA_SPACK_VERSION="1.2.2"
GYSELA_SPACK_PREFIX=${ALL_CCFRSCRATCH}/gysela-spack-${GYSELA_SPACK_VERSION}
export SPACK_DISABLE_LOCAL_CONFIG=true
export PYTHONDONTWRITEBYTECODE=True

if [ ! -d "${GYSELA_SPACK_PREFIX}" ]; then
    mkdir --parents "${GYSELA_SPACK_PREFIX}"
    chmod g+s "${GYSELA_SPACK_PREFIX}"
    git clone --branch v${GYSELA_SPACK_VERSION} --depth 1 https://github.com/spack/spack.git "${GYSELA_SPACK_PREFIX}"
fi

. ${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh

for arch in genoa mi250 mi300; do
    env="gyselalibxx-${arch}"
    env_file="${TOOLCHAIN_ROOT_DIRECTORY}/${arch}/gyselalibxx-spack-environment.yaml"

    echo "Preparing the Spack environment ${env}"

    spack env remove --yes-to-all "${env}"
    spack env create "${env}" "${env_file}"

    spack env activate "${env}"
    spack repo update
    spack concretize --quiet
    spack spec --install-status --namespaces
    spack install --jobs "${GYSELA_SPACK_JOBS}"
    spack env deactivate
done

find "${GYSELA_SPACK_PREFIX}" -user "${USER}" \
    -exec chgrp "${GYSELA_SPACK_GROUP}" {} + \
    -exec chmod g+rwX {} +
