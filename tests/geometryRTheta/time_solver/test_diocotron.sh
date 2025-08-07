#!/bin/bash
set -xe

if [ $# -ne 3 ]; then
    echo "Usage: $0 <GYSELALIBXX_SRCDIR> <GYSELALIBXX_EXEC> <PYTHON3_EXE>"
    exit 1
fi
GYSELALIBXX_SRCDIR="$1"
GYSELALIBXX_EXEC="$2"
PYTHON3_EXE="$3"

OUTDIR="${PWD}/${SIMULATION_NAME}"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
    rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
TESTDIR="${PWD}"

cd "${TMPDIR}"

# Create a parameter file with the default values. 
"${GYSELALIBXX_EXEC}" "--dump-config" "${PWD}/diocotron_params.yaml"

# Modify the default parameter files for a faster test.
MODIFY_PARAMETERS=$(cat <<EOF
import yaml

with open('diocotron_params.yaml') as f:
    data = yaml.safe_load(f)

data['SplineMesh']['r_ncells'] = 32
data['SplineMesh']['theta_ncells'] = 64
data['Time']['final_T'] = 40.
data['Output']['time_step_diag'] = 40

with open('diocotron_params.yaml', 'w') as f:
    yaml.dump(data, f)
EOF
)

python3 -c "$MODIFY_PARAMETERS"

# Launch the test with the modified parameter file. 
"${GYSELALIBXX_EXEC}" "${PWD}/diocotron_params.yaml"

export PYTHONPATH="${GYSELALIBXX_SRCDIR}/post-process/PythonScripts:${PYTHONPATH}"
"${PYTHON3_EXE}" -B "${GYSELALIBXX_SRCDIR}/tests/geometryRTheta/time_solver/growth_rate_test.py" . 
if [ ! -d ${OUTDIR} ]; then
    mkdir "${OUTDIR}"
fi
cp *.png "${OUTDIR}"