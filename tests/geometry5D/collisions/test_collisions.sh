#!/bin/bash
set -xe

# module load cray-python

if [ $# -ne 1 ]; then
    echo "Usage: $0 <INPUT_YAML_FILE>"
    exit 1
fi

# TESTCOLL_EXE="$(pwd)/../../../build.mi250.hipcc.adastra.spack/simulations/geometry5D/testcollisions/testcollisions"
TESTCOLL_EXE="$(pwd)/../../../build-v100/simulations/geometry5D/testcollisions/testcollisions"

GYSELA_IO_PATH="$(pwd)/gysela_io"

if [ ! -d ${GYSELA_IO_PATH} ]; then
    git clone git@gitlab.maisondelasimulation.fr:gysela-developpers/gysela_io.git
else
    (
        cd -- ${GYSELA_IO_PATH}
        git pull
    )
fi
export PYTHONPATH="${GYSELA_IO_PATH}:${PYTHONPATH}"

INPUT_YAML_FILE="${1}"

CASENAME=$(basename "${INPUT_YAML_FILE}" | tr '[:lower:]' '[:upper:]')
RESDIR="D_${CASENAME}"
if [ ! -d ${RESDIR} ]; then
    mkdir -p -- "${RESDIR}"
fi

INPUT_CXX_YAML_FILE='coll_ref.yml'
INIT_RST_FILE='GysX_rst_00000.h5'
SAVE_RST_FILE='GysX_rst_00001.h5'

cat <<EOF >${INPUT_CXX_YAML_FILE}
InputFileNames:
  read_restart : ${INIT_RST_FILE}
  write_restart : ${SAVE_RST_FILE}

Algorithm:
  deltat: 0.125
  nbiter: 360

Output:
  time_diag: 0.25
EOF

cp -- "${INPUT_YAML_FILE}" "${RESDIR}"
mv -- "${INPUT_CXX_YAML_FILE}" "${RESDIR}"
cd -- "${RESDIR}"
python3 "${GYSELA_IO_PATH}/initialisation/init_distribution.py" -i "${INPUT_YAML_FILE}" -o "${INIT_RST_FILE}"
if [ -f ${SAVE_RST_FILE} ]; then
  rm -- "${SAVE_RST_FILE}"
fi
${TESTCOLL_EXE} "${INPUT_CXX_YAML_FILE}"
python3 "${GYSELA_IO_PATH}/diagnostics/compare_slice.py" -i "${INIT_RST_FILE} ${SAVE_RST_FILE}" -s "tor1=0,tor2=0,tor3=0,species=0"

cd ..

echo "Results saved in ${RESDIR}."

sha1sum -- "${TESTCOLL_EXE}" "${RESDIR}/${INPUT_CXX_YAML_FILE}" "${RESDIR}/${INIT_RST_FILE}" "${RESDIR}/${SAVE_RST_FILE}"
