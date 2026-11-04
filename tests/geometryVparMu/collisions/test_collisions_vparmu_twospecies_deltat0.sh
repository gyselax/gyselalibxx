#!/bin/bash

set -eu
# set -o pipefail
# set -x

if [ ${#} -ne 1 ]; then
    echo "Usage: ${0} <COLL_EXEC>"
    exit 1
fi

SCRIPT_DIRECTORY="$(dirname -- "$(realpath -- "${BASH_SOURCE[0]:-${0}}")")"
TMP_DIRECTORY="$(mktemp -p "${SCRIPT_DIRECTORY}" -d run-XXXXXXXXXX)"
function finish {
    rm -rf -- "${TMP_DIRECTORY}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

COLL_EXEC="$(realpath -- "${1}")"

cd -- "${TMP_DIRECTORY}"

REFERENCE_DIRECTORY="${SCRIPT_DIRECTORY}/references"
if hostname | grep -q '^persee'; then
    REFERENCE_DIRECTORY="/data/gyselarunner/gyselalibxx_references/tests/geometryVparMu/collisions"
fi

LOG_FILE="${SCRIPT_DIRECTORY}/validation.log"
rm -f -- "${LOG_FILE}"
{
    hostname
    timedatectl
    ################################################################################
    # dt=0, checks the projection/reconstruction
    echo "1.delta_t_at_zero"
    mkdir 1.delta_t_at_zero && cd 1.delta_t_at_zero
    echo "Dump the default configuration coll2species_deltat0.yaml"
    "${COLL_EXEC}" --dump-config "coll2species_deltat0.yaml"
    sed "s/\(deltat:\) .*/\1 0./" "coll2species_deltat0.yaml" >"case.yaml"
    /usr/bin/time --verbose "${COLL_EXEC}" case.yaml 1>&2
    ../../validate.py coll_00000.h5 coll_00001.h5 fdistribu || true # Absorb errors for this one
    ../../validate.py coll_00001.h5 "${REFERENCE_DIRECTORY}/1.delta_t_at_zero/coll_00001.h5" fdistribu
    cd ..
    ################################################################################
    COLLISION_VPAR_MU_SIMULATIONS_PATH="${SCRIPT_DIRECTORY}/../../../simulations/geometryVparMu/collisions"
    # coll_1specie_maxwellian
    echo "2.coll_1specie_maxwellian"
    mkdir 2.coll_1specie_maxwellian && cd 2.coll_1specie_maxwellian
    /usr/bin/time --verbose "${COLL_EXEC}" "${COLLISION_VPAR_MU_SIMULATIONS_PATH}/coll_1specie_maxwellian.yaml" 1>&2
    ../../validate.py coll_00000.h5 "${REFERENCE_DIRECTORY}/2.coll_1specie_maxwellian/coll_00000.h5" fdistribu
    ../../validate.py coll_00010.h5 "${REFERENCE_DIRECTORY}/2.coll_1specie_maxwellian/coll_00010.h5" fdistribu
    cd ..
    ################################################################################
    # coll_2species_maxwellian_differentT
    echo "3.coll_2species_maxwellian_differentT"
    mkdir 3.coll_2species_maxwellian_differentT && cd 3.coll_2species_maxwellian_differentT
    /usr/bin/time --verbose "${COLL_EXEC}" "${COLLISION_VPAR_MU_SIMULATIONS_PATH}/coll_2species_maxwellian_differentT.yaml" 1>&2
    ../../validate.py coll_00000.h5 "${REFERENCE_DIRECTORY}/3.coll_2species_maxwellian_differentT/coll_00000.h5" fdistribu
    ../../validate.py coll_00010.h5 "${REFERENCE_DIRECTORY}/3.coll_2species_maxwellian_differentT/coll_00010.h5" fdistribu
    cd ..
    ################################################################################
    # coll_2species_maxwellian_differentV
    echo "4.coll_2species_maxwellian_differentV"
    mkdir 4.coll_2species_maxwellian_differentV && cd 4.coll_2species_maxwellian_differentV
    /usr/bin/time --verbose "${COLL_EXEC}" "${COLLISION_VPAR_MU_SIMULATIONS_PATH}/coll_2species_maxwellian_differentV.yaml" 1>&2
    ../../validate.py coll_00000.h5 "${REFERENCE_DIRECTORY}/4.coll_2species_maxwellian_differentV/coll_00000.h5" fdistribu
    ../../validate.py coll_00010.h5 "${REFERENCE_DIRECTORY}/4.coll_2species_maxwellian_differentV/coll_00010.h5" fdistribu
    cd ..
} | tee "${LOG_FILE}"
