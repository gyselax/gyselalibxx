#!/bin/bash
set -xe

if [ $# -ne 4 ]; then
    echo "Usage: $0 <GYSELALIBXX_SRCDIR> <GYSELALIBXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME>"
    exit 1
fi
GYSELALIBXX_SRCDIR="$1"
GYSELALIBXX_EXEC="$2"
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

"${GYSELALIBXX_EXEC}" "--dump-config" "${PWD}/diocotron_params.yaml"
"${GYSELALIBXX_EXEC}" "${PWD}/diocotron_params.yaml"

# # Theoretical values for Landau damping
# growthrate_theory=-0.153
# frequency_theory=1.4156

export PYTHONPATH="${GYSELALIBXX_SRCDIR}/post-process/PythonScripts:${PYTHONPATH}"
"${PYTHON3_EXE}" -B "${GYSELALIBXX_SRCDIR}/tests/geometryRTheta/time_solver/growth_rate_test.py" . #-g ${growthrate_theory} -f ${frequency_theory}
if [ ! -d ${OUTDIR} ]; then
    mkdir "${OUTDIR}"
fi
cp *.png "${OUTDIR}"