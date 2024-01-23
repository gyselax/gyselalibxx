#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "vortex_merger_equilibrium.hpp"



/**
 * @brief Initial condition for the vortex merger simulation.
 *
 * @tparam Mapping
 *      A Curvilinear2DToCartesian mapping function class.
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
     *      A mapping function from the logical domain
     *      to the physical domain.
     */
    VortexMergerDensitySolution(Mapping const& mapping) : m_mapping(mapping) {}

    ~VortexMergerDensitySolution() {};

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
            DSpanRP rho_init,
            DViewRP rho_eq,
            const double eps,
            const double sigma,
            const double x_star_1,
            const double y_star_1,
            const double x_star_2,
            const double y_star_2)
    {
        IDomainRP grid = ddc::get_domain<IDimR, IDimP>(rho_init);

        // Initialisation:
        for_each(grid, [&](IndexRP const irp) {
            const CoordRP coord_rp(ddc::coordinate(irp));
            const CoordXY coord_xy(m_mapping(coord_rp));
            const double x = ddc::get<RDimX>(coord_xy);
            const double y = ddc::get<RDimY>(coord_xy);

            rho_init(irp) = rho_eq(irp)
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
