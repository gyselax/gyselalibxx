if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

export OMP_PROC_BIND=false
export OMP_NUM_THREADS=4
