if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

. /data/gyselarunner/spack-0.20.0/share/spack/setup-env.sh

spack load gcc@11
spack env activate gyselalibxx-env-omp-cuda-ginkgo-1_8_0
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=8

