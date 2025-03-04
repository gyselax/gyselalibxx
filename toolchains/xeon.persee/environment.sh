if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

. /data/gyselarunner/spack-0.23.0/share/spack/setup-env.sh

spack load gcc@12
spack env activate gyselalibxx-env-omp
