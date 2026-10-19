#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    exit 1
fi

if [[ $# -lt 3 ]]; then
    echo "Usage: bash prepare.sh <env-base-name> <env-yaml-file> <env-conf-dir>" >&2
    exit 1
fi
GYSELA_ENV_BASE_NAME="$1"
GYSELA_ENV_YAML_FILE="$2"
GYSELA_ENV_CONF_DIR="$3"

set -eu

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

export PYTHONDONTWRITEBYTECODE=True

VERSION=$(cat ${TOOLCHAIN_ROOT_DIRECTORY}/../version.txt)

env="${GYSELA_ENV_BASE_NAME}--${VERSION}"
echo "Preparing the Spack environment ${env}"

spack env remove --yes-to-all "${env}"
spack env create "${env}"

env_dir="$(spack location --env ${env})"
cp "${GYSELA_ENV_YAML_FILE}" "${env_dir}/"
mkdir --parents "${env_dir}/gysela"
cp "${GYSELA_ENV_CONF_DIR}"/* "${env_dir}/gysela/"

spack env activate "${env}"
spack repo update
spack concretize --quiet
spack spec --install-status --long --namespaces --types
spack install --jobs "${GYSELA_SPACK_JOBS}"
spack env deactivate

find "${GYSELA_SPACK_PREFIX}" \( -type d -o -type f \) -user "${USER}" \
    -exec chgrp "${GYSELA_SPACK_GROUP}" {} + \
    -exec chmod g+rwX {} +
