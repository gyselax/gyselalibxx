#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

if [ -d "${HOME}/.spack" ]; then
    echo "It is recommended to remove any old spack installation $HOME/.spack before running this script."
    echo "$HOME/.spack detected, continue anyway? [YyNn]"
    read ignore_existing
    if [[ $ignore_existing == y* ]]; then
        echo "Ignoring `$HOME/.spack`"
    else
        exit 1
    fi
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Download spack
wget https://github.com/spack/spack/releases/download/v0.23.0/spack-0.23.0.tar.gz
tar -xvf spack-0.23.0.tar.gz
rm spack-0.23.0.tar.gz

# Bug fix for spack < 0.24 and gcc>=14
wget https://raw.githubusercontent.com/spack/spack/b369d8b2509794c4f46f62c81f25c247ca58418e/var/spack/repos/builtin/packages/py-netcdf4/package.py -O spack-0.23.0/var/spack/repos/builtin/packages/py-netcdf4/package.py

# Activate spack
. spack-0.23.0/share/spack/setup-env.sh

# Reduce the naming scheme of packages to avoid shebang issues.
# Increase the time out that is by default too short for some packages (like PDI)
spack config --scope site add 'config:connect_timeout:60'

# Restrict `blas` and `lapack` providers to avoid `amdlibflame` and `amdblis`, for some reason they might get selected and fail to compile.
spack config --scope site add 'packages:all:providers:blas:[openblas]'
spack config --scope site add 'packages:all:providers:lapack:[openblas]'

# Add PDI repository
git clone https://github.com/pdidev/spack.git spack-0.23.0/var/spack/repos/pdi
spack repo add --scope site spack-0.23.0/var/spack/repos/pdi

spack compiler find

AVAILABLE_COMPILERS=$(spack compilers | grep "gcc@1[1-9]" || true)

if [ -z "${AVAILABLE_COMPILERS}" ]
then
    echo "A gcc compiler with a version of at least 11 was not found. Installing compiler"
    spack install gcc@11
    spack load gcc@11
    spack compiler find
fi

spack env create gyselalibxx-env ${SCRIPT_DIR}/gyselalibxx-env-0.23.0.yaml
export SPACK_DEBUG_LOG_DIR=spack-logs
export SPACK_VERBOSE=1
export SPACK_DEBUG=1
spack --env gyselalibxx-env install --fail-fast --show-log-on-error --jobs 2


