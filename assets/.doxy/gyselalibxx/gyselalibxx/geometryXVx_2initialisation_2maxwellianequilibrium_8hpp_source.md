

# File maxwellianequilibrium.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**maxwellianequilibrium.hpp**](geometryXVx_2initialisation_2maxwellianequilibrium_8hpp.md)

[Go to the documentation of this file](geometryXVx_2initialisation_2maxwellianequilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iequilibrium.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

class MaxwellianEquilibrium : public IEquilibrium
{
    host_t<DFieldMemSp> m_density_eq;

    host_t<DFieldMemSp> m_temperature_eq;

    host_t<DFieldMemSp> m_mean_velocity_eq;

public:
    MaxwellianEquilibrium(
            host_t<DFieldMemSp> density_eq,
            host_t<DFieldMemSp> temperature_eq,
            host_t<DFieldMemSp> mean_velocity_eq);

    ~MaxwellianEquilibrium() override = default;

    static MaxwellianEquilibrium init_from_input(
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    DFieldSpVx operator()(DFieldSpVx allfequilibrium) const override;

    static void compute_maxwellian(
            DFieldVx const fMaxwellian,
            double const density,
            double const temperature,
            double const mean_velocity);

    host_t<DConstFieldSp> density_eq() const
    {
        return get_const_field(m_density_eq);
    }

    host_t<DConstFieldSp> temperature_eq() const
    {
        return get_const_field(m_temperature_eq);
    }

    host_t<DConstFieldSp> mean_velocity_eq() const
    {
        return get_const_field(m_mean_velocity_eq);
    }
};
```


