// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <typeinfo>

#include <sll/mapping/cartesian_to_circular.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>

#include "advection_domain.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "directional_tag.hpp"
#include "geometry.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "l_norm_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

/**
 * @brief Define a domain for the advection.
 *
 * The natural advection domain is the physical domain (AdvectionDomain),
 * where the studied equation is given.
 * However, not all the mappings used are analytically invertible and inverting
 * the Jacobian matrix of the mapping could be costly. That is why, we also
 * introduce a pseudo-Cartesian domain (AdvectionPseudoCartesianDomain).
 *
 * The method used to advect the characteristic feet depends on the
 * advection domain.
 *
 * More details can be found in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * @see BslAdvectionRTheta
 * @see IFootFinder
 */
template <class LogicalToPhysicalMapping>
class AdvectionDomain
{
    static_assert(is_analytical_mapping_v<LogicalToPhysicalMapping>);

private:
    using PhysicalToLogicalMapping = inverse_mapping_t<LogicalToPhysicalMapping>;

public:
    /**
     * @brief The first dimension in the advection domain.
     */
    using X_adv = typename PhysicalToLogicalMapping::cartesian_tag_x;
    /**
     * @brief The second dimension in the advection domain.
     */
    using Y_adv = typename PhysicalToLogicalMapping::cartesian_tag_y;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_adv = Coord<X_adv, Y_adv>;

private:
    LogicalToPhysicalMapping m_to_cartesian_mapping;
    PhysicalToLogicalMapping m_to_curvilinear_mapping;

public:
    /**
     * @brief Instantiate a AdvectionDomain advection domain.
     *
     * @param[in] to_physical_mapping
     *      The mapping from the logical domain to the physical domain.
     */
    AdvectionDomain(LogicalToPhysicalMapping const& to_physical_mapping)
        : m_to_cartesian_mapping(to_physical_mapping)
        , m_to_curvilinear_mapping(to_physical_mapping.get_inverse_mapping())
    {
    }

    /**
     * @brief Advect the characteristic feet.
     *
     * In the Backward Semi-Lagrangian method, the advection of a function
     * uses the conservation along the characteristic property. So, we firstly
     * compute the characteristic feet and then interpolate the function at these
     * characteristic feet.
     *
     * The function implemented here deals with the computation of the characteristic feet.
     * The IFootFinder class uses a time integration method to solve the characteristic
     * equation.
     * The BslAdvectionRTheta class calls advect_feet to compute the characteristic feet
     * and interpolate the function we want to advect.
     *
     * The advect_feet implemented here computes only
     *
     * - @f$  (\text{feet}_r, \text{feet}_\theta) =  \mathcal{F}^{-1} (\text{feet}_x, \text{feet}_y)  @f$
     *
     * with
     *      - @f$ \text{feet}_x = x_{ij} - dt A_x(x_{ij}, y_{ij}) @f$,
     *      - @f$ \text{feet}_y = y_{ij} - dt A_y(x_{ij}, y_{ij}) @f$,
     *
     * and
     *      - @f$ (x,y)_{ij} =   \mathcal{F}(r_{ij}, \theta_{ij})@f$, with @f$\{(r, \theta)_{ij}\}_{ij} @f$ the logical mesh points,
     *      - @f$ A @f$ the advection field in the advection domain,
     *      - @f$  \mathcal{F} @f$ the mapping from the logical domain to the advection domain.
     *
     * @param[in, out] feet_coords_rp
     *      The computed characteristic feet in the logical domain.
     *      On input: the points we want to advect.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field defined on the advection domain.
     * @param[in] dt
     *      The time step.
     */
    void advect_feet(
            host_t<FieldRTheta<CoordRTheta>> feet_coords_rp,
            host_t<DConstVectorFieldRTheta<X_adv, Y_adv>> advection_field,
            double dt) const
    {
        IdxRangeRTheta const idx_range_rp = get_idx_range<GridR, GridTheta>(feet_coords_rp);

        CoordXY_adv coord_center(m_to_cartesian_mapping(CoordRTheta(0, 0)));

        ddc::for_each(idx_range_rp, [&](IdxRTheta const irp) {
            CoordRTheta const coord_rp(feet_coords_rp(irp));
            CoordXY_adv const coord_xy = m_to_cartesian_mapping(coord_rp);

            CoordXY_adv const feet_xy = coord_xy - dt * advection_field(irp);

            if (norm_inf(feet_xy - coord_center) < 1e-15) {
                feet_coords_rp(irp) = CoordRTheta(0, 0);
            } else {
                feet_coords_rp(irp) = m_to_curvilinear_mapping(feet_xy);
                ddc::select<Theta>(feet_coords_rp(irp)) = ddcHelper::restrict_to_idx_range(
                        ddc::select<Theta>(feet_coords_rp(irp)),
                        IdxRangeTheta(idx_range_rp));
            }
        });
    }
};
