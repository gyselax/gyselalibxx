#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

current_folder=$(realpath $(dirname ${BASH_SOURCE[0]}))
gyselalibxx_folder=$(realpath ${current_folder}/../..)

cmake_prefixes="/opt/googletest:/opt/openmp/"

docker pull ghcr.io/gyselax/gyselalibxx_env:latest
docker run -v ${gyselalibxx_folder}:/src --workdir "/src" --user :$(id -g) -e CMAKE_PREFIX_PATH=${cmake_prefixes} -it ghcr.io/gyselax/gyselalibxx_env:latest

