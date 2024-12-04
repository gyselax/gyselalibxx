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
#include "math_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

/**
 * @brief Define a domain for the advection.
 *
 * The natural advection domain is the physical domain (AdvectionPhysicalDomain),
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
public:
    virtual ~AdvectionDomain() {};
};


/**
 * @brief Define the physical domain for the advection.
 *
 * The physical domain can only be used for analytically invertible mappings.
 *
 * Here the advection field in the physical domain and the advection domain is the same.
 * So the AdvectionPhysicalDomain::compute_advection_field return the same advection field.
 *
 * The method to advect is here given by:
 * - compute @f$ (x,y)_{ij} = \mathcal{F} (r,\theta)_{ij} @f$,
 * - advect
 *      - @f$ \text{feet}_{x, ij} = x_{ij} - dt A_x(x_{ij}, y_{ij}) @f$,
 *      - @f$ \text{feet}_{y, ij} = y_{ij} - dt A_y(x_{ij}, y_{ij}) @f$,
 * - compute @f$ (\text{feet}_{r, ij}, \text{feet}_{\theta, ij})
 * = \mathcal{F}^{-1} (\text{feet}_{x, ij}, \text{feet}_{y, ij}) @f$
 * as @f$\mathcal{F} @f$ easily invertible.
 *
 *
 */
template <class LogicalToPhysicalMapping, class PhysicalToLogicalMapping>
class AdvectionPhysicalDomain : public AdvectionDomain<LogicalToPhysicalMapping>
{
public:
    /**
     * @brief The first dimension in the advection domain.
     */
    using X_adv = X;
    /**
     * @brief The second dimension in the advection domain.
     */
    using Y_adv = Y;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_adv = Coord<X_adv, Y_adv>;

private:
    LogicalToPhysicalMapping const& m_to_cartesian_mapping;
    PhysicalToLogicalMapping const& m_to_curvilinear_mapping;

public:
    /**
     * @brief Instantiate a AdvectionPhysicalDomain advection domain.
     *
     * @param[in] to_physical_mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] to_logical_mapping
     *      The mapping from the physical domain to the logical domain.
     */
    AdvectionPhysicalDomain(
            LogicalToPhysicalMapping const& to_physical_mapping,
            PhysicalToLogicalMapping const& to_logical_mapping)
        : m_to_cartesian_mapping(to_physical_mapping)
        , m_to_curvilinear_mapping(to_logical_mapping) {};
    ~AdvectionPhysicalDomain() {};

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
        using namespace ddc;

        IdxRangeRTheta const idx_range_rp = get_idx_range<GridR, GridTheta>(feet_coords_rp);
        CoordXY coord_center(m_to_cartesian_mapping(CoordRTheta(0, 0)));

        ddc::for_each(idx_range_rp, [&](IdxRTheta const irp) {
            CoordRTheta const coord_rp(feet_coords_rp(irp));
            CoordXY const coord_xy = m_to_cartesian_mapping(coord_rp);

            CoordXY const feet_xy = coord_xy - dt * advection_field(irp);

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


/**
 * @brief Define the pseudo-Cartesian domain for the advection.
 *
 * The pseudo-Cartesian domain is recommended for not analytically
 * invertible mappings.
 *
 * We introduce a circular mapping @f$ \mathcal{G}@f$
 * from the logical domain to the pseudo-Cartesian domain.
 * The circular mapping is analytically invertible.
 *
 * The advection field is defined in the physical domain, so we need to
 * define it in the pseudo-Cartesian domain before advecting
 * (AdvectionPseudoCartesianDomain::compute_advection_field).
 * To do so, we use the Jacobian matrix @f$ J_{\mathcal{F}}J_{\mathcal{G}}^{-1} @f$.
 * This matrix is invertible :
 *
 * @f$ A_{\text{cart}} = (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1} A @f$.
 *
 * The main difficulty is the center point. So for @f$ r \in [0, \varepsilon] @f$, we
 * linearize the Jacobian matrix:
 *
 * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1} (r, \theta)
 * = \left(1 - \frac{r}{\varepsilon} \right) (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1} (0, \theta)
 * +  \frac{r}{\varepsilon}   (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1} (\varepsilon, \theta)
 * @f$.
 *
 * The @f$  (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1} (0, \theta) @f$ matrices are implemented
 * in the Curvilinear2DToCartesian::to_pseudo_cartesian_jacobian_center_matrix child classes.
 *
 *
 * The method to advect is here given by:
 * - compute @f$ (x_{\text{cart}},y_{\text{cart}})_{ij} = \mathcal{G} (r,\theta)_{ij} @f$,
 * - advect
 *      - @f$ \text{feet}_{x_{\text{cart}}, ij} = x_{\text{cart}, ij} - dt A_{\text{cart}, x}(x_{\text{cart}, ij},  y_{\text{cart}, ij}) @f$,
 *      - @f$ \text{feet}_{y_{\text{cart}}, ij} = y_{\text{cart}, ij} - dt A_{\text{cart}, y}(x_{\text{cart}, ij},  y_{\text{cart}, ij}) @f$,
 * - compute @f$ (\text{feet}_{r, ij}, \text{feet}_{\theta, ij})
 * = \mathcal{G}^{-1} (\text{feet}_{x_{\text{cart}}, ij}, \text{feet}_{y_{\text{cart}}, ij}) @f$.
 *
 *
 *
 * More details can be found in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * @see DiscreteToCartesian
 */
template <class LogicalToPhysicalMapping>
class AdvectionPseudoCartesianDomain : public AdvectionDomain<LogicalToPhysicalMapping>
{
public:
    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix2x2 = std::array<std::array<double, 2>, 2>;

    /**
     * @brief The first dimension in the advection domain.
     */
    using X_adv = DimX_pC;
    /**
     * @brief The second dimension in the advection domain.
     */
    using Y_adv = DimY_pC;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_adv = Coord<X_adv, Y_adv>;


private:
    LogicalToPhysicalMapping const& m_to_cartesian_mapping;
    double m_epsilon;

public:
    /**
     * @brief Instantiate an AdvectionPseudoCartesianDomain advection domain.
     *
     * @param[in] to_physical_mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] epsilon
     *      @f$ \varepsilon @f$ parameter used for the linearization of the
     *      advection field around the central point.
     */
    AdvectionPseudoCartesianDomain(
            LogicalToPhysicalMapping const& to_physical_mapping,
            double epsilon = 1e-12)
        : m_to_cartesian_mapping(to_physical_mapping)
        , m_epsilon(epsilon) {};
    ~AdvectionPseudoCartesianDomain() {};

    /**
     * @brief Advect the characteristic feet.
     *
     * In the Backward Semi-Lagrangian method, the advection of a function
     * uses the conservation along the characteristic property. So, we firstly
     * compute the characteristic feet and then interpolate the function at this
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
     * - @f$  (\text{feet}_r, \text{feet}_\theta) =  \mathcal{G}^{-1} (\text{feet}_x, \text{feet}_y)  @f$
     *
     * with
     *      - @f$ \text{feet}_x = x_{ij} - dt A_x(x_{ij}, y_{ij}) @f$,
     *      - @f$ \text{feet}_y = y_{ij} - dt A_y(x_{ij}, y_{ij}) @f$,
     *
     * and
     *      - @f$ (x,y)_{ij} =   \mathcal{G}(r_{ij}, \theta_{ij})@f$, with @f$\{(r, \theta)_{ij}\}_{ij} @f$ the logical mesh points,
     *      - @f$ A @f$ the advection field in the advection domain,
     *      - @f$  \mathcal{G} @f$ the mapping from the logical domain to the advection domain.
     *
     * @param[in, out] feet_coords_rp
     *      The computed characteristic feet in the logical domain.
     *      On input: the points we want to advect.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field defined in the advection domain.
     * @param[in] dt
     *      The time step.
     */
    void advect_feet(
            host_t<FieldRTheta<CoordRTheta>> feet_coords_rp,
            host_t<DConstVectorFieldRTheta<X_adv, Y_adv>> const& advection_field,
            double const dt) const
    {
        static_assert(
                !std::is_same_v<LogicalToPhysicalMapping, CircularToCartesian<R, Theta, X, Y>>);
        IdxRangeRTheta const idx_range_rp = get_idx_range(advection_field);

        CircularToCartesian<R, Theta, X_adv, Y_adv> const pseudo_Cartesian_mapping;
        CoordXY_adv const center_xy_pseudo_cart
                = CoordXY_adv(pseudo_Cartesian_mapping(CoordRTheta(0., 0.)));

        ddc::for_each(idx_range_rp, [&](IdxRTheta const irp) {
            CoordRTheta const coord_rp(feet_coords_rp(irp));
            CoordXY_adv const coord_xy_pseudo_cart = pseudo_Cartesian_mapping(coord_rp);
            CoordXY_adv const feet_xy_pseudo_cart
                    = coord_xy_pseudo_cart - dt * advection_field(irp);

            if (norm_inf(feet_xy_pseudo_cart - center_xy_pseudo_cart) < 1e-15) {
                feet_coords_rp(irp) = CoordRTheta(0, 0);
            } else {
                CartesianToCircular<X_adv, Y_adv, R, Theta> const inv_pseudo_Cartesian_mapping;
                feet_coords_rp(irp) = inv_pseudo_Cartesian_mapping(feet_xy_pseudo_cart);
                ddc::select<Theta>(feet_coords_rp(irp)) = ddcHelper::restrict_to_idx_range(
                        ddc::select<Theta>(feet_coords_rp(irp)),
                        IdxRangeTheta(idx_range_rp));
            }
        });
    }
};
