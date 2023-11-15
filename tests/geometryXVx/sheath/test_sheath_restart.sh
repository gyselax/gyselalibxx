#!/bin/bash
set -xe

if [ $# -ne 3 ]
then
    echo "Usage: $0 <VOICEXX_SRCDIR> <VOICEXX_EXEC> <SIMULATION_NAME>"
    exit 1
fi
VOICEXX_SRCDIR="$1"
VOICEXX_EXEC="$2"
SIMULATION_NAME="$3"

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

"${VOICEXX_EXEC}" "--dump-config" "${PWD}/sheath.yaml"
sed -i 's/^  x_size: .*/  x_size: 16/' sheath.yaml
sed -i 's/^  vx_size: .*/  vx_size: 16/' sheath.yaml
sed -i 's/^  nbiter: .*/  nbiter: 10/' sheath.yaml
sed -i 's/^  deltat: .*/  deltat: 0.125/' sheath.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.25/' sheath.yaml

"${VOICEXX_EXEC}" "${PWD}/sheath.yaml"

cd "${TMPDIR}"
cp "${RSTDIR}/VOICEXX_initstate.h5" .
cp "${RSTDIR}/VOICEXX_00003.h5" .
cp "${RSTDIR}/sheath.yaml" sheath_restart.yaml

"${VOICEXX_EXEC}" --iter-restart 3 "${PWD}/sheath_restart.yaml"
sed -i 's/^  nbiter: .*/  nbiter: 2/' sheath_restart.yaml
sed -i 's/^  time_diag: .*/  time_diag: 0.5/' sheath_restart.yaml

h5ls -d ${PWD}/VOICEXX_00006.h5/time_saved ${RSTDIR}/VOICEXX_00005.h5/time_saved
command="h5diff ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 time_saved"
eval $command
if [ $? -ne 0 ]; then
    exit 1
fi

command="h5diff ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 electrostatic_potential"
eval $command
if [ $? -ne 0 ]; then
    exit 1
fi

command="h5diff ${PWD}/VOICEXX_00005.h5 ${RSTDIR}/VOICEXX_00005.h5 fdistribu"
eval $command
if [ $? -ne 0 ]; then
    exit 1
fi