# Libraries in Gyselalib++

## Documentation of the libraries

- DDC (discrete domain computation): [https://ddc.mdls.fr/index.html](https://ddc.mdls.fr/index.html).
- Eigen: [https://eigen.tuxfamily.org/index.php?title=Main_Page#Documentation](https://eigen.tuxfamily.org/index.php?title=Main_Page#Documentation).
- Gingko (linear algrebra): [https://github.com/ginkgo-project/ginkgo/wiki/Tutorial:-Building-a-Poisson-Solver](https://github.com/ginkgo-project/ginkgo/wiki/Tutorial:-Building-a-Poisson-Solver).
- GoogleTest (unit test): [https://google.github.io/googletest/](https://google.github.io/googletest/).
- Kokkos (parallelization GPU): [https://kokkos.org/](https://kokkos.org/). 
- Koliop (collision operator): [https://gitlab.com/cines/code.gysela/libkoliop](https://gitlab.com/cines/code.gysela/libkoliop). 
- PDI (input/output): [https://pdi.dev/main/](https://pdi.dev/main/). 



## Additional examples of use 
### Save and load data with PDI
PDI is a library giving an interface to read and write input/output for simulations. 
Here is a simple example of use. More information can be found on the documentation 
[https://pdi.dev/main/](https://pdi.dev/main/). 

In the `main.cpp` file where we want to save the data, 

```cpp
#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>
#include <paraconf.h>
#include <pdi.h>

#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"

int main(int argc, char** argv)
{
    PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Read from a .ymal parameter file
    double const dt(PCpp_double(conf_gyselalibxx, ".Time.delta_t")); // example for a double.
    int const time_step(PCpp_int(conf_gyselalibxx, ".Time.time_step")); // example for an integer.

    // Write in a .h5 output file
    ddc::expose_to_pdi("delta_t", dt); // for double/integer.
    expose_mesh_to_pdi("x_coords", mesh_x); // for meshes (see src/io)

    ddc::PdiEvent("initialization") // save data on this event 
        .with("data_1", data_1) // 1D array (see pdi_out.yml.hpp)
        .with("data_2", data_2) // 2D array (see pdi_out.yml.hpp)
        .with("data_3", data_3); // 3D array (see pdi_out.yml.hpp)
    // The data can be double/integer/array on host

    for (int iter(0); iter < 10; iter++){
        ddc::PdiEvent("iter") // save data on this event 
            .with("data_1", data_1);
    }
}
```

In the same folder as the `main.cpp` file, add `pdi_out.yml.hpp` file with the details of the data types:
```cpp
// SPDX-License-Identifier: MIT
#pragma once

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  delta_t: double 
  time_step: int

  x_coords_extents: {type: array,  subtype: int64, size: 1 }
  x_coords:
    type: array
    subtype: double
    size: [  '$x_coords_extents[0]' ]

data:
  data_1_extents: {type: array,  subtype: int64, size: 1 }
  data_1:
    type: array
    subtype: double
    size: [  '$data_1_extents[0]' ]

  data_2_extents: {type: array,  subtype: int64, size: 2 }
  data_2:
    type: array
    subtype: double
    size: [  '$data_2_extents[0]', '$data_2_extents[1]' ]

  data_3_extents: {type: array,  subtype: int64, size: 3 }
  data_3:
    type: array
    subtype: double
    size: [  '$data_3_extents[0]', '$data_3_extents[1]', '$data_3_extents[2]' ]


plugins:
  set_value:
    on_init:
      - share:
        - iter_saved: 0
    on_data:
      iteration:
        - set:
          - iter_saved: '${iter} % ${2}' 
    on_finalize:
      - release: [iter_saved]

  decl_hdf5:
    - file: 'output/VOICEXX_initstate.h5'   // name of the saving file.
      on_event: [initialization]            
      collision_policy: replace_and_warn    // if the file exist, replace the file and warn. 
      write: [data_1, data_2, data_3]       // list of data we save on the event. 

    - file: 'output/VOICEXX_${iter:05}.h5'   // name of the saving file.
      on_event: [iteration, last_iteration] 
      when: '${iter} % ${2} = 0'             // save the data at this condition
      collision_policy: replace_and_warn
      write: [data_1]
  #trace: ~
)PDI_CFG";

```

Once the simulation has run, we treat the data with the python files in the `post-process` folder.
We use `gysdata` to load the data. Each geometry has its one data structure. Be careful that the type you define
in the `pdi_out.yml.hpp` are also present in the `data_structure_XXX.yaml` you use. 
In the python file, get the saved data using `DiskStore` as follows. 
```py 
from pathlib import Path
from gysdata import DiskStore

path_data_structure = Path('data_structure_XXX.yaml')
ds = DiskStore(folder, data_structure=path_data_structure)

dt = ds["delta_t"]
time_step = ds["time_step"]

data_1 = np.array(ds["data_1"])
x_coords = np.array(ds["data_1"].coords("x").values)
```

The last line corresponds to the case where `data_1` where defined on `x_coords` in the 
`data_structure_XXX.yaml` file. 
```yaml
# one field
data_1:
  # the list of HDF5 datasets where this field can be found (this can also be a single dict)
  path:
    # file regex, in the regex we capture values that can be used in filename_coord (?P<NAME>CAPTURE)
    file: 'VOICEXX_\d+.h5'
    # dataset regex
    dataset: 'data_1'
  dimensions:
  # ordered set of dimension_name -> how to get the coordinates of this dimension
  - &xd
    x:       { global_coord: [ &initstate "VOICEXX_initstate.h5", x_coords ] }
```


