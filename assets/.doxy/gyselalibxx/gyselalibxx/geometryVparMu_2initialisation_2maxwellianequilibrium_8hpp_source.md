

# File maxwellianequilibrium.hpp

[**File List**](files.md) **>** [**geometryVparMu**](dir_9a2f28dc8f538ee0f4428810facf29b8.md) **>** [**initialisation**](dir_99d29839093a8e7b0be0d596be7efa54.md) **>** [**maxwellianequilibrium.hpp**](geometryVparMu_2initialisation_2maxwellianequilibrium_8hpp.md)

[Go to the documentation of this file](geometryVparMu_2initialisation_2maxwellianequilibrium_8hpp.md)


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
    // mass of all kinetic species
    host_t<DFieldMemSp> m_mass;

    // equilibrium density of all kinetic species
    host_t<DFieldMemSp> m_density_eq;

    // equilibrium temperature of all kinetic species
    host_t<DFieldMemSp> m_temperature_eq;

    // equilibrium mean velocity of all kinetic species
    host_t<DFieldMemSp> m_mean_velocity_eq;

private:
    // magnetic field
    double m_magnetic_field;

public:
    MaxwellianEquilibrium(
            host_t<DFieldMemSp> mass,
            host_t<DFieldMemSp> density_eq,
            host_t<DFieldMemSp> temperature_eq,
            host_t<DFieldMemSp> mean_velocity_eq,
            double magnetic_field);

    ~MaxwellianEquilibrium() override = default;

    static MaxwellianEquilibrium init_from_input(
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    DFieldSpVparMu operator()(DFieldSpVparMu allfequilibrium) const override;


    static void compute_maxwellian(
            DFieldVparMu const fMaxwellian,
            double const mass,
            double const density,
            double const temperature,
            double const mean_velocity,
            double const magnetic_field);

    host_t<DConstFieldSp> mass() const
    {
        return get_const_field(m_mass);
    }

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


