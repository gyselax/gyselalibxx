#!/bin/bash
set -xe

if [ $# -lt 4 ] || [ $# -gt 6 ]
then
    echo "Usage: $0 <GYSELALIBXX_SRCDIR> <GYSELALIBXX_EXEC> <PYTHON3_EXE> <SIMULATION_NAME> [<RELATIVE_TOLERANCE> <ABSOLUTE_TOLERANCE>]"
    exit 1
fi
GYSELALIBXX_SRCDIR="$1"
GYSELALIBXX_EXEC="$2"
PYTHON3_EXE="$3"
SIMULATION_NAME="$4"
if [ -n "$5" ]
then
  RELATIVE_TOLERANCE="$5"
fi
if [ -n "$6" ]
then
  ABSOLUTE_TOLERANCE="$6"
fi

TMPDIR="$(mktemp -p "${PWD}" -d run-XXXXXXXXXX)"

#function finish {
#  rm -rf "${TMPDIR}"
#}
#trap finish EXIT QUIT ABRT KILL SEGV TERM STOP

cd ${TMPDIR}

git clone git@gitlab.maisondelasimulation.fr:gysela-developpers/gysela_io.git

${PYTHON3_EXE} ${TMPDIR}/gysela_io/initialisation/init_distribution.py -i ${TMPDIR}/gysela_io/initialisation/input_examples/input_params_ref4D.yaml -o ${TMPDIR}/GyselaX_twospecies16x32x128x64_00000.h5

cd ${TMPDIR}

"${GYSELALIBXX_EXEC}" --dump-config "${TMPDIR}/config.yml"

"${GYSELALIBXX_EXEC}" "${TMPDIR}/config.yml"

${PYTHON3_EXE} ${GYSELALIBXX_SRCDIR}/post-process/PythonScripts/compare_hdf5_results.py "${TMPDIR}/GyselaX_twospecies16x32x128x64_00000.h5" "${TMPDIR}/GyselaX_twospecies16x32x128x64_00001.h5" fdistribu -R ${RELATIVE_TOLERANCE} -A ${ABSOLUTE_TOLERANCE}

