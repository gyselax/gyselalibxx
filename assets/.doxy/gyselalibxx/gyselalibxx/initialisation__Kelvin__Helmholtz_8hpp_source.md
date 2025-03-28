

# File initialisation\_Kelvin\_Helmholtz.hpp

[**File List**](files.md) **>** [**geometryXY**](dir_8727f3a3f911772a0d72e99b040a604a.md) **>** [**initialisation**](dir_bbd6cbd2460631a86c53df5e1524a2e1.md) **>** [**initialisation\_Kelvin\_Helmholtz.hpp**](initialisation__Kelvin__Helmholtz_8hpp.md)

[Go to the documentation of this file](initialisation__Kelvin__Helmholtz_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"


class KelvinHelmholtzInstabilityInitialisation
{
    double const m_epsilon;
    double const m_mode_k;

public:
    KelvinHelmholtzInstabilityInitialisation(double const epsilon, double const mode_k)
        : m_epsilon(epsilon)
        , m_mode_k(mode_k) {};

    ~KelvinHelmholtzInstabilityInitialisation() = default;

    void operator()(DFieldXY allfdistribu, DFieldXY allfdistribu_equilibrium)
    {
        IdxRangeXY const idx_range = get_idx_range(allfdistribu);
        double const epsilon_proxy = m_epsilon;
        double const mode_k_proxy = m_mode_k;

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const i_xy) {
                    CoordXY const coord_xy(ddc::coordinate(i_xy));
                    double const x = CoordX(coord_xy);
                    double const y = CoordY(coord_xy);
                    allfdistribu_equilibrium(i_xy) = Kokkos::sin(y);
                    allfdistribu(i_xy) = allfdistribu_equilibrium(i_xy)
                                         + epsilon_proxy * Kokkos::cos(mode_k_proxy * x);
                });
    };
};
```


