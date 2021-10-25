#!/bin/bash
set -xe

if [ $# -ne 3 ]
then
    echo "Usage: $0 <VOICEXX_SRCDIR> <VOICEXX_EXEC> <PYTHON3_EXE>"
    exit 1
fi
VOICEXX_SRCDIR="$1"
VOICEXX_EXEC="$2"
PYTHON3_EXE="$3"

OUTDIR="${PWD}"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
  rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
TESTDIR="${PWD}"

cd "${TMPDIR}"

"${VOICEXX_EXEC}" "${TESTDIR}/params.yaml"
export PYTHONPATH="${VOICEXX_SRCDIR}/post-process/PythonScripts"
"${PYTHON3_EXE}" -B "${TESTDIR}/check_Landau.py"
cp *.png "${OUTDIR}"
