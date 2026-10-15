#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

GYSELA_SPACK_GROUP="ksw"
GYSELA_SPACK_JOBS=128
GYSELA_SPACK_VERSION="1.2.2"
GYSELA_SPACK_BASE=${ALL_CCFRSCRATCH}
