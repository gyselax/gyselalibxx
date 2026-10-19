#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    exit 1
fi

# check that GYSELA_SPACK_PREFIX is set

if [ ! -d "${GYSELA_SPACK_PREFIX}" ]; then
    mkdir --parents "${GYSELA_SPACK_PREFIX}"
    chmod g+s "${GYSELA_SPACK_PREFIX}"
    git clone --branch v${GYSELA_SPACK_VERSION} --depth 1 https://github.com/spack/spack.git "${GYSELA_SPACK_PREFIX}"
fi
