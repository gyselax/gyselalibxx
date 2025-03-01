// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "vortex_merger_equilibrium.hpp"



/**
 * @brief Initial condition for the vortex merger simulation.
 *
 * @tparam Mapping
 *      A class describing a mapping from curvilinear coordinates to Cartesian coordinates.
 */
template <class Mapping>
class VortexMergerDensitySolution
{
private:
    Mapping const& m_mapping;


public:
    /**
     * @brief Instantiate a VortexMergerDensitySolution.
     *
     * @param[in] mapping
     *      A mapping function from the logical index range
     *      to the physical index range.
     */
    explicit VortexMergerDensitySolution(Mapping const& mapping) : m_mapping(mapping) {}

    /**
     * @brief Set an initial condition for the vortex merger simulation.
     *
     * The initial condition is given by
     * @f$ \rho(0, x, y) = \rho_{eq}(x,y) + \varepsilon
     *      \left(
     *          \exp\left[-\frac{(x - x_1^*)^2 + (y - y_1^*)^2}{2\sigma^2}\right]
     *       +  \exp\left[-\frac{(x - x_2^*)^2 + (y - y_2^*)^2}{2\sigma^2}\right]
     *      \right)@f$
     *
     * with @f$ \rho_{eq}@f$ given by VortexMergerEquilibria::set_equilibrium.
     *
     *
     * The initial condition is also given in Edoardo Zoni's article
     * (https://doi.org/10.1016/j.jcp.2019.108889).
     *
     * @param[out] rho_init
     *      The initial condition of the density.
     * @param[in] rho_eq
     *      The equilibrium density solution computed thanks to
     *      VortexMergerEquilibria::set_equilibrium.
     * @param[in] eps
     *      The @f$ \varepsilon @f$ amplitude of the perturbation.
     * @param[in] sigma
     *      The @f$ \sigma @f$ of the Gaussian functions.
     * @param[in] x_star_1
     *      The @f$ x_1^*@f$ parameter.
     * @param[in] y_star_1
     *      The @f$ y_1^*@f$ parameter.
     * @param[in] x_star_2
     *      The @f$ x_2^*@f$ parameter.
     * @param[in] y_star_2
     *      The @f$ y_2^*@f$ parameter.
     *
     */
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
