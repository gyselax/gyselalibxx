#!/bin/bash

[[ -z ${1+x} ]] && echo -e "Gyselalibxx preparation shortcut script:\n    ./prepare.gyselalibxx.sh <toolchain/machine [mi250.cce.adastra.spack|...]> <do prepare [default: TRUE|FALSE]> <force a fresh build [default: TRUE|FALSE]> <extra CMake flags>" && exit 1

PREPARE_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

set -eu

PREPARE_TOOLCHAIN_NAME="${1}"

if [[ ! -z ${2+x} ]]; then
    PREPARE_DONT_SKIP="${2}"
else
    PREPARE_DONT_SKIP="TRUE"
fi

if [[ ! -z ${3+x} ]]; then
    PREPARE_FORCE_FRESH_BUILD="${3}"
else
    PREPARE_FORCE_FRESH_BUILD="TRUE"
fi

# A toolchain script should NOT expect we source/execute them from the toolchain
# directory.
PREPARE_TOOLCHAIN_PATH="${PREPARE_ROOT_DIRECTORY}/toolchains/${PREPARE_TOOLCHAIN_NAME}"

if [[ "${PREPARE_DONT_SKIP}" != "TRUE" ]]; then
    echo "[PREPARE] Skipping the environment preparation (presumably already prepared)."
else
    # Anything that does not match "TRUE" disable the feature.
    echo "[PREPARE] Preparing the environment. Preparation log in 'log.prepare.${PREPARE_TOOLCHAIN_NAME}'."
    # Maybe a noop.
    "${PREPARE_TOOLCHAIN_PATH}/prepare.sh" >"log.prepare.${PREPARE_TOOLCHAIN_NAME}" 2>&1 || (echo "The environment preparation failed! Check the log." && exit 1)
fi

PREPARE_BUILD_NAME="build.${PREPARE_TOOLCHAIN_NAME}"
PREPARE_BUILD_DIRECTORY="${PREPARE_ROOT_DIRECTORY}/${PREPARE_BUILD_NAME}"
PREPARE_TOOLCHAIN_FILE="${PREPARE_TOOLCHAIN_PATH}/toolchain.cmake"

if [[ "${PREPARE_FORCE_FRESH_BUILD}" != "TRUE" ]]; then
    echo "[PREPARE] Keeping the old build directory."
else
    # Anything that does not machine "TRUE" disable the feature.
    echo "[PREPARE] Removing old build directory..."
    rm -rf -- "${PREPARE_BUILD_DIRECTORY}"
fi

PREPARE_BUILD_GENERATOR="Unix Makefiles"

echo "[PREPARE] Sourcing environment file: '${PREPARE_TOOLCHAIN_PATH}/environment.sh'."

source -- "${PREPARE_TOOLCHAIN_PATH}/environment.sh"

if which ninja >/dev/null 2>&1; then
    PREPARE_BUILD_GENERATOR="Ninja"
fi

echo "[PREPARE] Using build system: '${PREPARE_BUILD_GENERATOR}'."

PREPARE_BUILD_CMAKE_FLAG="${@:4}"

cmake \
    -B "${PREPARE_BUILD_DIRECTORY}" \
    -S . \
    -G "${PREPARE_BUILD_GENERATOR}" \
    -DCMAKE_TOOLCHAIN_FILE="${PREPARE_TOOLCHAIN_FILE}" \
    "${PREPARE_BUILD_CMAKE_FLAG}"

echo "Building in '${PREPARE_BUILD_NAME}'. Build log in 'log.build.${PREPARE_TOOLCHAIN_NAME}'."
cmake --build "${PREPARE_BUILD_DIRECTORY}" >"log.build.${PREPARE_TOOLCHAIN_NAME}" 2>&1 || (echo "The build failed! Check the log." && exit 1)

echo "The build completed. You can use the binaries in the following directory:"
echo "    \"${PREPARE_BUILD_DIRECTORY}\""
echo "Do not forget to source the environment file like so:"
echo "    $ source -- \"${PREPARE_TOOLCHAIN_PATH}/environment.sh\""
