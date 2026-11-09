#!/bin/bash
set -xe

if [ $# -ne 5 ]
then
    echo "Usage: $0 <GYSELALIBXX_SRCDIR> <GYSELALIBXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME> <INPUT_LANDAU>"
    exit 1
fi
GYSELALIBXX_SRCDIR="$1"
GYSELALIBXX_EXEC="$2"
PYTHON3_EXE="$3"
SIMULATION_NAME="$4"
INPUT_LANDAU="$5"

OUTDIR="${PWD}/${SIMULATION_NAME}"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
  rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
TESTDIR="${PWD}"

cd "${TMPDIR}"
cp ${INPUT_LANDAU} . 
"${GYSELALIBXX_EXEC}" "${INPUT_LANDAU}"



