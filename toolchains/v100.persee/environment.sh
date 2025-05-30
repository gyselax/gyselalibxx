if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

if [ -f "$(which spack)" ]
then
  spack env deactivate
else
  . /data/gyselarunner/spack-0.23.0/share/spack/setup-env.sh
fi

spack load gcc@12
spack env activate gyselalibxx-env-omp-cuda
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=8

