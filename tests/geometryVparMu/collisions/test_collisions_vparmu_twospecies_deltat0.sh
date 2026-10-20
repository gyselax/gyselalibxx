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

# Generate the input YAML file by default
INPUT_YAML_FILE="${PWD}/coll2species_deltat0.yaml"
echo "Create the input YAML file by default ${INPUT_YAML_FILE}"
"${COLL_EXEC}" --dump-config "${INPUT_YAML_FILE}"
# Change deltat to 0. in the input YAML file.
NEW_DELTAT=0.
sed -i "s/\(deltat:\) .*/\1 ${NEW_DELTAT}/" "${INPUT_YAML_FILE}"
"${COLL_EXEC}" "${INPUT_YAML_FILE}"
../validate.py coll_00000.h5 coll_00001.h5 fdistribu |& tee ../references/xyz.log
echo "Success"
