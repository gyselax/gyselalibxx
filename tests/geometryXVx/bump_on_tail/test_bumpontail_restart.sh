#!/bin/bash
set -xe

if [ $# -lt 4 ] || [ $# -gt 6 ]
then
    echo "Usage: $0 <GYSELALIBXX_SRCDIR> <GYSELALIBXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME> [<RELATIVE_RESTART_TOLERANCE> <ABSOLUTE_RESTART_TOLERANCE>]"
    exit 1
fi
GYSELALIBXX_SRCDIR="$1"
GYSELALIBXX_EXEC="$2"
PYTHON3_EXE="$3"
SIMULATION_NAME="$4"
if [ -n "$5" ]
then
  RELATIVE_RESTART_TOLERANCE="$5"
fi
if [ -n "$6" ]
then
  ABSOLUTE_RESTART_TOLERANCE="$6"
fi

OUTDIR="${PWD}/${SIMULATION_NAME}"

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"
function finish {
  rm -rf "${TMPDIR}"
}
trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd "$(dirname "$0")"
TESTDIR="${PWD}"

cd "${TMPDIR}"

RSTDIR="${TMPDIR}/RST"
mkdir "${RSTDIR}"
cd "${RSTDIR}"

"${GYSELALIBXX_EXEC}" "--dump-config" "${PWD}/bumpontail.yaml"
sed -i 's/^  x_size: .*/  x_size: 16/' bumpontail.yaml
sed -i 's/^  vx_size: .*/  vx_size: 16/' bumpontail.yaml
sed -i 's/^  nbiter: .*/  nbiter: 10/' bumpontail.yaml
sed -i 's/^  deltat: .*/  deltat: 0.125/' bumpontail.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.25/' bumpontail.yaml

"${GYSELALIBXX_EXEC}" "${PWD}/bumpontail.yaml"

cd "${TMPDIR}"
cp "${RSTDIR}/GYSELALIBXX_initstate.h5" .
cp "${RSTDIR}/GYSELALIBXX_00003.h5" .
cp "${RSTDIR}/bumpontail.yaml" bumpontail_restart.yaml

"${GYSELALIBXX_EXEC}" --iter-restart 3 "${PWD}/bumpontail_restart.yaml"
sed -i 's/^  nbiter: .*/  nbiter: 2/' bumpontail_restart.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.5/' bumpontail_restart.yaml

h5ls -d ${PWD}/GYSELALIBXX_00005.h5/time_saved ${RSTDIR}/GYSELALIBXX_00005.h5/time_saved
${PYTHON3_EXE} ${GYSELALIBXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/GYSELALIBXX_00005.h5 ${RSTDIR}/GYSELALIBXX_00005.h5 time_saved -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi

${PYTHON3_EXE} ${GYSELALIBXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/GYSELALIBXX_00005.h5 ${RSTDIR}/GYSELALIBXX_00005.h5 electrostatic_potential -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi

${PYTHON3_EXE} ${GYSELALIBXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/GYSELALIBXX_00005.h5 ${RSTDIR}/GYSELALIBXX_00005.h5 fdistribu -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi
