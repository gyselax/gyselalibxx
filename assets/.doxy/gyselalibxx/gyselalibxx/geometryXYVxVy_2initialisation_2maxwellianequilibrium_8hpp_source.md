

# File maxwellianequilibrium.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**initialisation**](dir_51031f497920158ed20948cdaeaff0bc.md) **>** [**maxwellianequilibrium.hpp**](geometryXYVxVy_2initialisation_2maxwellianequilibrium_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2initialisation_2maxwellianequilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry_xyvxvy.hpp"
#include "iequilibrium.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

class MaxwellianEquilibrium : public IEquilibrium
{
    // equilibrium density of all kinetic species
    host_t<FieldMemSp<Real>> m_density_eq;

    // equilibrium temperature of all kinetic species
    host_t<FieldMemSp<Real>> m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    host_t<FieldMemSp<Real>> m_mean_velocity_eq;

public:
    MaxwellianEquilibrium(
            host_t<FieldMemSp<Real>> density_eq,
            host_t<FieldMemSp<Real>> temperature_eq,
            host_t<FieldMemSp<Real>> mean_velocity_eq);

    ~MaxwellianEquilibrium() override = default;

    static MaxwellianEquilibrium init_from_input(
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    DFieldSpVxVy operator()(DFieldSpVxVy allfequilibrium) const override;


    static void compute_maxwellian(
            DFieldVxVy const fMaxwellian,
            Real const density,
            Real const temperature,
            Real const mean_velocity);

    host_t<ConstFieldSp<Real>> density_eq() const
    {
        return get_const_field(m_density_eq);
    }

    host_t<ConstFieldSp<Real>> temperature_eq() const
    {
        return get_const_field(m_temperature_eq);
    }

    host_t<ConstFieldSp<Real>> mean_velocity_eq() const
    {
        return get_const_field(m_mean_velocity_eq);
    }
};
```


