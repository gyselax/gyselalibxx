// SPDX-License-Identifier: MIT

#pragma once
#include <geometry.hpp>

#include "ddc_aliases.hpp"


/**
 * @brief Initialize the allfdistribu function. 
 *  
 * Set the allfdistribu at 
 * @f$ f_{eq}(x, y) = \sin(y) @f$
 * 
 * and 
 * @f$ f_0(x, y) = f_{eq}(x, y) + \varepsilon \cos(kx) @f$, 
 * 
 * with @f$ \varepsilon @f$ an amplitude of perturbation and 
 * @f$ k @f$ mode equal to @f$ 2\pi @f$ divided by the length of the 
 * index range on @f$ x @f$. 
 */
class KelvinHelmholtzInstabilityInitialization
{
    double const m_epsilon;
    double const m_mode_k;

public:
    /**
     * @brief Instantiate the initializer. 
     * 
     * @param epsilon @f$ \varepsilon @f$, the amplitude of perturbation. 
     * @param mode_k @f$ k @f$, the perturbation mode. 
     */
    KelvinHelmholtzInstabilityInitialization(double const epsilon, double const mode_k)
        : m_epsilon(epsilon)
        , m_mode_k(mode_k) {};

    ~KelvinHelmholtzInstabilityInitialization() = default;

    /**
     * @brief Initialise @f$ f_{eq}@f$ and @f$ f @f$.
     * 
     * @param allfdistribu Field refering to the @f$ f @f$ function.
     * @param allfdistribu_equilibrium Field refering to the @f$ f_{eq} @f$ function.
     */
    void operator()(DFieldXY allfdistribu, DFieldXY allfdistribu_equilibrium)
    {
        IdxRangeXY const idx_range = get_idx_range(allfdistribu);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const i_xy) {
                    CoordXY const coord_xy(ddc::coordinate(i_xy));
                    double const x = CoordX(coord_xy);
                    double const y = CoordY(coord_xy);
                    allfdistribu_equilibrium(i_xy) = Kokkos::sin(y);
                    allfdistribu(i_xy) = allfdistribu_equilibrium(i_xy)
                                         + m_epsilon * Kokkos::cos(m_mode_k * x);
                });
    };
};
