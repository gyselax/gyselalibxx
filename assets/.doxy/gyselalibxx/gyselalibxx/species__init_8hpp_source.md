

# File species\_init.hpp

[**File List**](files.md) **>** [**speciesinfo**](dir_661be8452a62f1b4720eb6eb57123ae7.md) **>** [**species\_init.hpp**](species__init_8hpp.md)

[Go to the documentation of this file](species__init_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <paraconf.h>

#include "ddc_helper.hpp"
#include "paraconfpp.hpp"
#include "pdi_helper.hpp"
#include "species_info.hpp"

void init_all_species(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_gyselalibxx,
        int nb_kinspecies,
        int nb_fluidspecies);

IdxRangeSp init_species(PC_tree_t conf_gyselalibxx);

void init_species_withfluid(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_gyselalibxx);

IdxRangeSp init_kinetic_species();
```


