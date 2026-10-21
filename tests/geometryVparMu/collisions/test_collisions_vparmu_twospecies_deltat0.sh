#!/bin/bash

set -eu
# set -x

if [ $# -ne 1 ]; then
    echo "Usage: $0 <COLL_EXEC>"
    exit 1
fi
COLL_EXEC="${1}"
TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
    rm -rf -- "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd -- "${TMPDIR}"
INPUT_YAML_FILE="${PWD}/coll2species_deltat0.yaml"
echo "Dump the default configuration: '${INPUT_YAML_FILE}'"
"${COLL_EXEC}" --dump-config "${INPUT_YAML_FILE}"

LOG_FILE="${PWD}/../references/mi250x.hipcc.adastra.spack.log"
rm -f -- "${LOG_FILE}"
################################################################################
# dt=0, checks the projection/reconstruction
mkdir 1.delta_t_at_zero && cd 1.delta_t_at_zero
sed "s/\(deltat:\) .*/\1 0./" "${INPUT_YAML_FILE}" >case.yaml
"${COLL_EXEC}" case.yaml
../../validate.py coll_00000.h5 coll_00001.h5 fdistribu |& tee -a "${LOG_FILE}"
cd ..
