# Collisions test

<!-- Script to automise :
 - the creation of the initial restart file with the python script `init_distribution.py`
 - the creation of the input YAML file required as input of the C++ simulation `testcollision`-->

<!-- For instance, the following command -->
`./testcollisions.sh input_params_ion_1x1x1x128x64.yaml`

<!-- create the folder `D_INPUT_PARAMS_ION_1X1X1X128X64.YAML` containing:
  - `GysX_rst_00000.h5` : output of the python script `init_distribution.py`
  - `GysX_rst_00001.h5` : output of the C++ collision executable
  - `coll_ref.yml` : input for C++ collision executable automatically created by the bash script `testcollision.sh`
  - `diff_f_vpar_mu_itor1eq0_itor2eq0_itor3eq0_ispeq0.png` : output figure to compare the results between `GysX_rst_00000.h5` and `GysX_rst_00001.h5`-->

