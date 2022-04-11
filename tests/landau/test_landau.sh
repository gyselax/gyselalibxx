#!/bin/bash
set -xe

if [ $# -ne 4 ]
then
    echo "Usage: $0 <VOICEXX_SRCDIR> <VOICEXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME>"
    exit 1
fi
VOICEXX_SRCDIR="$1"
VOICEXX_EXEC="$2"
PYTHON3_EXE="$3"
SIMULATION_NAME="$4"

OUTDIR="${PWD}/${SIMULATION_NAME}"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
  rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
TESTDIR="${PWD}"

cd "${TMPDIR}"

"${VOICEXX_EXEC}" "--dump-config" "${TESTDIR}/landau.yaml"
"${VOICEXX_EXEC}" "${TESTDIR}/landau.yaml"
export PYTHONPATH="${VOICEXX_SRCDIR}/post-process/PythonScripts"
"${PYTHON3_EXE}" -B "${TESTDIR}/check_Landau.py"
if [ ! -d ${OUTDIR} ]
then
    mkdir "${OUTDIR}"
fi
cp *.png "${OUTDIR}"

