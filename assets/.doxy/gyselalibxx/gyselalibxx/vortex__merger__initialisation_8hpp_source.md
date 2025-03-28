

# File vortex\_merger\_initialisation.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**initialisation**](dir_1b70d60e6147eeeade38d183e3e9d318.md) **>** [**vortex\_merger\_initialisation.hpp**](vortex__merger__initialisation_8hpp.md)

[Go to the documentation of this file](vortex__merger__initialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "vortex_merger_equilibrium.hpp"



template <class Mapping>
class VortexMergerDensitySolution
{
private:
    Mapping const& m_mapping;


public:
    explicit VortexMergerDensitySolution(Mapping const& mapping) : m_mapping(mapping) {}

    void set_initialisation(
            host_t<DFieldRTheta> rho_init,
            host_t<DConstFieldRTheta> rho_eq,
            const double eps,
            const double sigma,
            const double x_star_1,
            const double y_star_1,
            const double x_star_2,
            const double y_star_2)
    {
        IdxRangeRTheta grid = get_idx_range<GridR, GridTheta>(rho_init);

        // Initialisation:
        ddc::for_each(grid, [&](IdxRTheta const irtheta) {
            const CoordRTheta coord_rtheta(ddc::coordinate(irtheta));
            const CoordXY coord_xy(m_mapping(coord_rtheta));
            const double x = ddc::get<X>(coord_xy);
            const double y = ddc::get<Y>(coord_xy);

            rho_init(irtheta) = rho_eq(irtheta)
                                + eps
                                          * (std::exp(
                                                     -((x - x_star_1) * (x - x_star_1)
                                                       + (y - y_star_1) * (y - y_star_1))
                                                     / (2 * sigma * sigma))
                                             + std::exp(
                                                     -((x - x_star_2) * (x - x_star_2)
                                                       + (y - y_star_2) * (y - y_star_2))
                                                     / (2 * sigma * sigma)));
        });
    }
};
```


