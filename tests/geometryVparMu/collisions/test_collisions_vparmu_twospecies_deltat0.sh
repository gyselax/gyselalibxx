#!/bin/bash
set -xe

if [ $# -ne 3 ]; then
    echo "Usage: $0 <COLL_SRCDIR> <COLL_EXEC> <PYTHON3_EXE>"
    exit 1
fi
COLL_SRCDIR="$1"
COLL_EXEC="$2"
PYTHON3_EXE="$3"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
    rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
cd "${TMPDIR}"

# --> Generate the input YAML file by default
INPUT_YAML_FILE="${PWD}/coll2species_deltat0.yaml"
echo "Create the input YAML file by default" $INPUT_YAML_FILE
"${COLL_EXEC}" "--dump-config" "${INPUT_YAML_FILE}"

# --> Change deltat to 0. in the input YAML file
new_deltat=0.
sed -i "s/\(deltat:\) .*/\1 $new_deltat/" "${INPUT_YAML_FILE}"

# --> Run the code
"${COLL_EXEC}" "${INPUT_YAML_FILE}"

# --> Compute the maximum error of the distribution function between initialisation and final time
ABSOLUTE_ERROR_TOLERANCE=1.e-1
FILE1="coll_00000.h5"
FILE2="coll_00001.h5"
DATASET="fdistribu"
h5dump -d $DATASET $FILE1 > temp1.txt
h5dump -d $DATASET $FILE2 > temp2.txt
MAX_ERROR=$(paste temp1.txt temp2.txt | awk '{diff=($1>$2) ? $1-$2 : $2-$1;print diff}' | awk 'BEGIN {max=0} {if ($1>max) max=$1} END {print max}')
rm temp1.txt temp2.txt

if awk "BEGIN {exit !($MAX_ERROR >$ABSOLUTE_ERROR_TOLERANCE)}";then
  echo "Error: Maximum error = $MAX_ERROR exceeds threshold $ABSOLUTE_ERROR_TOLERANCE"
  exit 1
fi

