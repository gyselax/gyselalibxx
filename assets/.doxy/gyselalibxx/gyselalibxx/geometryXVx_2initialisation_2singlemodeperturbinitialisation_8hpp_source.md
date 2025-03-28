

# File singlemodeperturbinitialisation.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**singlemodeperturbinitialisation.hpp**](geometryXVx_2initialisation_2singlemodeperturbinitialisation_8hpp.md)

[Go to the documentation of this file](geometryXVx_2initialisation_2singlemodeperturbinitialisation_8hpp.md)


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
    DConstFieldSpVx m_fequilibrium;

    host_t<IFieldMemSp> m_init_perturb_mode;

    host_t<DFieldMemSp> m_init_perturb_amplitude;

public:
    void perturbation_initialisation(
            DFieldX perturbation,
            int const perturb_mode,
            double const perturb_amplitude) const;

    SingleModePerturbInitialisation(
            DConstFieldSpVx fequilibrium,
            host_t<IFieldMemSp> init_perturb_mode,
            host_t<DFieldMemSp> init_perturb_amplitude);

    ~SingleModePerturbInitialisation() override = default;

    static SingleModePerturbInitialisation init_from_input(
            DConstFieldSpVx allfequilibrium,
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const override;
};
```


