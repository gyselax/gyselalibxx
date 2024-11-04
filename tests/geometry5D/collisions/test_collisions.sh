#!/bin/bash
set -xe

# Collisions test
#
#Script to automise :
# - the creation of the initial restart file with the python script `init_distribution.py`
# - the creation of the input YAML file required as input of the C++ simulation `testcollision`
#
#For instance, the following command
#`./testcollisions.sh ${GYSELALIBXX_SRC}/build/simulations/geometry5D/testcollisions input_params_twospecies_geom5D.yaml`
#
#creates the folder `D_INPUT_PARAMS_TWOSPECIES_GEOM5D.YAML` containing:
#  - `GysX_rst_00000.h5` : output of the python script `init_distribution.py`
#  - `GysX_rst_00001.h5` : output of the C++ collision executable
#  - `coll_ref.yml` : input for C++ collision executable automatically created by the bash script `testcollision.sh`
#  - `diff_f_vpar_mu_itor1eq0_itor2eq0_itor3eq0_ispeq0.png` : output figure to compare the results between `GysX_rst_00000.h5` and `GysX_rst_00001.h5`


# module load cray-python

if [ $# -ne 2 ]; then
    echo "Usage: $0 <TESTCOLLISION_EXE> <INPUT_YAML_FILE>"
    exit 1
fi

TESTCOLL_EXE=$(realpath $1)

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

INPUT_YAML_FILE="$2"

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
python3 "${GYSELA_IO_PATH}/initialisation/init_distribution.py"  -i "${INPUT_YAML_FILE}" -o "${INIT_RST_FILE}"
if [ -f ${SAVE_RST_FILE} ]; then
  rm -- "${SAVE_RST_FILE}"
fi
${TESTCOLL_EXE} "${INPUT_CXX_YAML_FILE}"
python3 "${GYSELA_IO_PATH}/diagnostics/compare_slice.py" -i "${INIT_RST_FILE} ${SAVE_RST_FILE}" -s "tor1=0,tor2=0,tor3=0,species=0"

cd ..

echo "Results saved in ${RESDIR}."

sha1sum -- "${TESTCOLL_EXE}" "${RESDIR}/${INPUT_CXX_YAML_FILE}" "${RESDIR}/${INIT_RST_FILE}" "${RESDIR}/${SAVE_RST_FILE}"
