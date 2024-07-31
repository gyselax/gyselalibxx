// SPDX-License-Identifier: MIT

#pragma once

#include <sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>

#include <directional_tag.hpp>
#include <geometry.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "advection_domain.hpp"
#include "i_interpolator_2d_rp.hpp"
#include "iadvectionrp.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"



/**
 * @brief Define an advection operator on 2D @f$(r, \theta)@f$ domain.
 *
 * The advection operator uses a semi-Lagrangian method. The method is based on
 * the property that the solution is constant along the characteristics.
 *
 * For the following equation:
 * @f$\partial_t f(t,x) + V(t, x) \cdot \nabla_x f(t,x) = 0,  @f$
 *
 * we write the characteristics:
 * @f$ \partial_t X(t; s, x) = V(t, X(t; s, x)), \qquad \text{ with } X(s; s, x) = x. @f$
 *
 * Then the property gives us:
 * @f$ f(t, x) = f(0, X(t; 0, x)), \quad \forall t. @f$
 *
 *
 * So the first step of the advection operator is to compute the characteristic feet @f$ X(t; t+dt, x_i) @f$
 * for each mesh point @f$ x_i @f$.
 *
 * For the second step, we interpolate the function at the characteristic feet computed, and obtain the
 * function at the next time step: @f$ f(t + dt, x) = f(t, X(t + dt; t, x))@f$.
 *
 *
 * Different time integration methods are implemented to solve the characteristic equation.
 * They are defined in the IFootFinder class.
 *
 * The feet can be advected in different domains.
 * Theses domains are defined in the AdvectionDomain class.
 *
 * The interpolation of the function is always done in the logical domain.
 *
 *
 *
 * @see IFootFinder
 * @see AdvectionDomain
 *
 */
template <class FootFinder, class Mapping>
class BslAdvectionRP : public IAdvectionRP
{
private:
    PreallocatableSplineInterpolatorRP<ddc::NullExtrapolationRule> const& m_interpolator;

    FootFinder const& m_find_feet;

    Mapping const& m_mapping;


public:
    /**
     * @brief Instantiate an advection operator.
     *
     * @param [in] function_interpolator
     *       The polar interpolator to interpolate the function once the
     *      characteristic computed.
     * @param[in] foot_finder
     *      An IFootFinder which computes the characteristic feet.
     * @param[in] mapping
     *      The mapping function from the logical domain to the physical
     *      domain. 
     *
     * @tparam IFootFinder
     *      A child class of IFootFinder.
     */
    BslAdvectionRP(
            PreallocatableSplineInterpolatorRP<ddc::NullExtrapolationRule> const&
                    function_interpolator,
            FootFinder const& foot_finder,
            Mapping const& mapping)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
        , m_mapping(mapping)
    {
    }

    ~BslAdvectionRP() override = default;


    /**
     * @brief Allocate a ChunkSpan to the advected function.
     *
     * @param [in, out] allfdistribu
     *      A ChunkSpan containing the values of the function we want to advect.
     * @param [in] advection_field_xy
     *      A VectorDViewRP containing the values of the advection field
     *      on the physical domain axes.
     * @param [in] dt
     *      A time step used.
     *
     * @return A ChunkSpan to allfdistribu advected on the time step given.
     */
    DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimX, RDimY> advection_field_xy,
            double dt) const
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRP> const interpolator_ptr = m_interpolator.preallocate();

        // Initialise the feet
        FieldRP<CoordRP> feet_rp(advection_field_xy.domain());
        ddc::for_each(advection_field_xy.domain(), [&](IndexRP const irp) {
            feet_rp(irp) = ddc::coordinate(irp);
        });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(feet_rp.span_view(), advection_field_xy, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, feet_rp.span_cview());

        return allfdistribu;
    }


    /**
     * @brief Allocate a ChunkSpan to the advected function.
     *
     * @param [in, out] allfdistribu
     *      A ChunkSpan containing the values of the function we want to advect.
     * @param [in] advection_field_rp
     *      A VectorDViewRP containing the values of the advection field
     *      on the logical domain axis.
     * @param [in] advection_field_xy_center
     *      A CoordXY containing the value of the advection field on the 
     *      physical domain axis at the O-point. 
     * @param [in] dt
     *      A time step used.
     *
     * @return A ChunkSpan to allfdistribu advected on the time step given.
     */
    DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimR, RDimP> advection_field_rp,
            CoordXY const& advection_field_xy_center,
            double dt) const
    {
        IDomainRP grid(allfdistribu.domain<IDimR, IDimP>());

        const int npoints_p = IDomainP(grid).size();
        IDomainRP const grid_without_Opoint(grid.remove_first(IVectRP(1, 0)));
        IDomainRP const Opoint_grid(grid.take_first(IVectRP(1, npoints_p)));


        // Convert advection field on RP to advection field on XY
        VectorDFieldRP<RDimX, RDimY> advection_field_xy(grid);

        ddc::for_each(grid_without_Opoint, [&](IndexRP const irp) {
            CoordRP const coord_rp(ddc::coordinate(irp));

            std::array<std::array<double, 2>, 2> J; // Jacobian matrix
            m_mapping.jacobian_matrix(coord_rp, J);
            std::array<std::array<double, 2>, 2> G; // Metric tensor
            m_mapping.metric_tensor(coord_rp, G);

            ddcHelper::get<RDimX>(advection_field_xy)(irp)
                    = ddcHelper::get<RDimR>(advection_field_rp)(irp) * J[1][1] / std::sqrt(G[1][1])
                      + ddcHelper::get<RDimP>(advection_field_rp)(irp) * -J[1][0]
                                / std::sqrt(G[0][0]);
            ddcHelper::get<RDimY>(advection_field_xy)(irp)
                    = ddcHelper::get<RDimR>(advection_field_rp)(irp) * -J[0][1] / std::sqrt(G[1][1])
                      + ddcHelper::get<RDimP>(advection_field_rp)(irp) * J[0][0]
                                / std::sqrt(G[0][0]);
        });

        ddc::for_each(Opoint_grid, [&](IndexRP const irp) {
            ddcHelper::get<RDimX>(advection_field_xy)(irp) = CoordX(advection_field_xy_center);
            ddcHelper::get<RDimY>(advection_field_xy)(irp) = CoordY(advection_field_xy_center);
        });

        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRP> const interpolator_ptr = m_interpolator.preallocate();

        // Initialise the feet
        FieldRP<CoordRP> feet_rp(grid);
        ddc::for_each(grid, [&](IndexRP const irp) { feet_rp(irp) = ddc::coordinate(irp); });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(feet_rp.span_view(), advection_field_xy, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, feet_rp.span_cview());

        return allfdistribu;
    }
};
