
if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

module purge
module load gcc/12.2.0

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

. ${TOOLCHAIN_ROOT_DIRECTORY}/spack_0.22.2/share/spack/setup-env.sh

spack load gcc@12
spack env activate ${HOME}/gyselalib-env
