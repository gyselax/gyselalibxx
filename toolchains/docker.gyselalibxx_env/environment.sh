
if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

current_folder=$(realpath $(dirname ${BASH_SOURCE[0]}))
gyselalibxx_folder=$(realpath ${current_folder}/../..)

docker run -v ${gyselalibxx_folder}:/src -it --entrypoint bash ghcr.io/gyselax/gyselalibxx_env -c "cd /src && bash"

