

# File bumpontailequilibrium.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**bumpontailequilibrium.hpp**](bumpontailequilibrium_8hpp.md)

[Go to the documentation of this file](bumpontailequilibrium_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iequilibrium.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

class BumpontailEquilibrium : public IEquilibrium
{
    host_t<DFieldMemSp> m_epsilon_bot;

    host_t<DFieldMemSp> m_temperature_bot;

    host_t<DFieldMemSp> m_mean_velocity_bot;

public:
    void compute_twomaxwellian(
            DFieldVx fMaxwellian,
            double epsilon_bot,
            double temperature_bot,
            double mean_velocity_bot) const;
    BumpontailEquilibrium(
            host_t<DFieldMemSp> epsilon_bot,
            host_t<DFieldMemSp> temperature_bot,
            host_t<DFieldMemSp> mean_velocity_bot);

    ~BumpontailEquilibrium() override = default;

    static BumpontailEquilibrium init_from_input(
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    DFieldSpVx operator()(DFieldSpVx allfequilibrium) const override;

    host_t<DConstFieldSp> epsilon_bot() const
    {
        return get_const_field(m_epsilon_bot);
    }

    host_t<DConstFieldSp> temperature_bot() const
    {
        return get_const_field(m_temperature_bot);
    }

    host_t<DConstFieldSp> mean_velocity_bot() const
    {
        return get_const_field(m_mean_velocity_bot);
    }
};
```


