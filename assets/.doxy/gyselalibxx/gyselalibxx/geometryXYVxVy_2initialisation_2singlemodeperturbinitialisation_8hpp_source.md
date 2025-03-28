

# File singlemodeperturbinitialisation.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**initialisation**](dir_51031f497920158ed20948cdaeaff0bc.md) **>** [**singlemodeperturbinitialisation.hpp**](geometryXYVxVy_2initialisation_2singlemodeperturbinitialisation_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2initialisation_2singlemodeperturbinitialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iinitialisation.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

class SingleModePerturbInitialisation : public IInitialisation
{
    DConstFieldSpVxVy m_fequilibrium;

    host_t<IFieldMemSp> m_init_perturb_mode;

    host_t<DFieldMemSp> m_init_perturb_amplitude;

public:
    void perturbation_initialisation(
            DFieldXY perturbation,
            int const perturb_mode,
            double const perturb_amplitude) const;

    SingleModePerturbInitialisation(
            DConstFieldSpVxVy fequilibrium,
            host_t<IFieldMemSp> init_perturb_mode,
            host_t<DFieldMemSp> init_perturb_amplitude);

    ~SingleModePerturbInitialisation() override = default;

    DFieldSpXYVxVy operator()(DFieldSpXYVxVy allfdistribu) const override;

    static SingleModePerturbInitialisation init_from_input(
            DConstFieldSpVxVy allfequilibrium,
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);
};
```


