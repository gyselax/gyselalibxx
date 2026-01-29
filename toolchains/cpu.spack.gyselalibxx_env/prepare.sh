#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

# set paths for spack to not use home ~/.spack folder
CURRENT_DIR=$(pwd)
export SPACK_PATH=${CURRENT_DIR}/spack-1.1.0/
export SPACK_USER_CONFIG_PATH="${SPACK_PATH}/user_config"
export SPACK_SYSTEM_CONFIG_PATH="${SPACK_PATH}/sys_config"
export SPACK_USER_CACHE_PATH="${SPACK_PATH}/user_cache"

set -e

if [ -d "${SPACK_PATH}" ]; then
    echo "It is recommended to remove any old spack installation ${SPACK_PATH} before running this script."
    echo "${SPACK_PATH} detected, continue anyway? [YyNn]"
    read ignore_existing
    if [[ $ignore_existing == y* ]]; then
        echo "Using existing ${SPACK_PATH}"
    else
        exit 1
    fi
else
    # Download spack
    wget https://github.com/spack/spack/releases/download/v1.1.0/spack-1.1.0.tar.gz
    tar -xf spack-1.1.0.tar.gz
    rm spack-1.1.0.tar.gz
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Activate spack
. ${SPACK_PATH}/share/spack/setup-env.sh

# Reduce the naming scheme of packages to avoid shebang issues.
# Increase the time out that is by default too short for some packages (like PDI)
spack config --scope site add 'config:connect_timeout:60'

spack compiler find

if [[ "$OSTYPE" == "darwin"* ]]; then
  COMPILER='apple-clang@14:'
else
  AVAILABLE_COMPILERS=$(spack compilers | grep "gcc@1[1-9]" || true)

  if [ -z "${AVAILABLE_COMPILERS}" ]
  then
      echo "A gcc compiler with a version of at least 11 was not found. Installing compiler"
      spack install gcc@11
      spack load gcc@11
      spack compiler find
  fi

  COMPILER='gcc@11:'
fi

spack env create gyselalibxx-env ${SCRIPT_DIR}/gyselalibxx-env-1.1.0.yaml
spack --env gyselalibxx-env install --jobs 2
spack env activate -p gyselalibxx-env
PYTHON_EXECUTABLE=$(which python3)
spack env deactivate

cat >${SCRIPT_DIR}/environment.sh <<EOL
if [ "\${BASH_SOURCE[0]}" -ef "\$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

export SPACK_PATH=${CURRENT_DIR}/spack-1.1.0/
export SPACK_USER_CONFIG_PATH="\${SPACK_PATH}/user_config"
export SPACK_SYSTEM_CONFIG_PATH="\${SPACK_PATH}/sys_config"
export SPACK_USER_CACHE_PATH="\${SPACK_PATH}/user_cache"
# The hdf5 package is injecting the environment view `lib` path to `LD_LIBRARY_PATH`
# which causes spurious segfaults for system executables, we manually remove it.
LD_LIBRARY_PATH_TMP="$LD_LIBRARY_PATH"
. \${SPACK_PATH}/share/spack/setup-env.sh
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH_TMP"
unset LD_LIBRARY_PATH_TMP
spack env activate -p gyselalibxx-env
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
EOL
