#!/bin/bash
set -xe

if [ $# -lt 4 ] || [ $# -gt 6 ]
then
    echo "Usage: $0 <VOICEXX_SRCDIR> <VOICEXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME> [<RELATIVE_RESTART_TOLERANCE> <ABSOLUTE_RESTART_TOLERANCE>]"
    exit 1
fi
VOICEXX_SRCDIR="$1"
VOICEXX_EXEC="$2"
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

${PYTHON3_EXE} ${VOICEXX_SRCDIR}/pre-process/PythonScripts/geometryXVx/suggested_points_refinement.py ${RSTDIR}/grids.h5 --edge-domains 0.0 50.0 --ncells 16 --name grid_x --periodic --xmin 0.0 --xmax 50.0
${PYTHON3_EXE} ${VOICEXX_SRCDIR}/pre-process/PythonScripts/geometryXVx/suggested_points_refinement.py ${RSTDIR}/grids.h5 --edge-domains -6.0 6.0 --ncells 16 --name grid_vx

"${VOICEXX_EXEC}" "--dump-config" "${PWD}/sheath.yaml"
sed -i 's/^  x_ncells: .*/  x_ncells: 16/' sheath.yaml
sed -i 's/^  vx_ncells: .*/  vx_ncells: 16/' sheath.yaml
sed -i 's/^  nbiter: .*/  nbiter: 10/' sheath.yaml
sed -i 's/^  deltat: .*/  deltat: 0.125/' sheath.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.25/' sheath.yaml

"${VOICEXX_EXEC}" "${PWD}/sheath.yaml"

cd "${TMPDIR}"
cp "${RSTDIR}/VOICEXX_initstate.h5" .
cp "${RSTDIR}/VOICEXX_00003.h5" .
cp "${RSTDIR}/sheath.yaml" sheath_restart.yaml
cp "${RSTDIR}/grids.h5" .

"${VOICEXX_EXEC}" --iter-restart 3 "${PWD}/sheath_restart.yaml"
sed -i 's/^  nbiter: .*/  nbiter: 2/' sheath_restart.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.5/' sheath_restart.yaml

h5ls -d ${PWD}/VOICEXX_00005.h5/time_saved ${RSTDIR}/VOICEXX_00005.h5/time_saved
${PYTHON3_EXE} ${VOICEXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 time_saved -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi

${PYTHON3_EXE} ${VOICEXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 electrostatic_potential -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi

${PYTHON3_EXE} ${VOICEXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 fdistribu -R ${RELATIVE_RESTART_TOLERANCE} -A ${ABSOLUTE_RESTART_TOLERANCE}
if [ $? -ne 0 ]; then
    exit 1
fi
