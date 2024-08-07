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
#include "ddc_aliases.hpp"
#include "i_interpolator_2d_rp.hpp"
#include "iadvectionrp.hpp"
#include "spline_foot_finder.hpp"
#include "spline_interpolator_2d_rp.hpp"



/**
 * @brief Define an advection operator on 2D @f$(r, \theta)@f$ index range.
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
 * The feet can be advected in different index ranges.
 * Theses index ranges are defined in the AdvectionDomain class.
 *
 * The interpolation of the function is always done in the logical index range.
 *
 *
 *
 * @see IFootFinder
 * @see AdvectionDomain
 *
 */
template <class FootFinder, class Mapping>
class BslAdvectionRTheta : public IAdvectionRTheta
{
private:
    PreallocatableSplineInterpolatorRTheta<ddc::NullExtrapolationRule> const& m_interpolator;

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
     *      The mapping function from the logical index range to the physical
     *      index range. 
     *
     * @tparam IFootFinder
     *      A child class of IFootFinder.
     */
    BslAdvectionRTheta(
            PreallocatableSplineInterpolatorRTheta<ddc::NullExtrapolationRule> const&
                    function_interpolator,
            FootFinder const& foot_finder,
            Mapping const& mapping)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
        , m_mapping(mapping)
    {
    }

    ~BslAdvectionRTheta() override = default;


    /**
     * @brief Allocate a Field of the advected function.
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_xy
     *      A DConstVectorFieldRTheta containing the values of the advection field
     *      on the physical index range axes.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<X, Y> advection_field_xy,
            double dt) const
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRTheta> const interpolator_ptr = m_interpolator.preallocate();

        // Initialise the feet
        FieldMemRTheta<CoordRTheta> feet_rp(get_idx_range(advection_field_xy));
        ddc::for_each(get_idx_range(advection_field_xy), [&](IdxRTheta const irp) {
            feet_rp(irp) = ddc::coordinate(irp);
        });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(get_field(feet_rp), advection_field_xy, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, get_const_field(feet_rp));

        return allfdistribu;
    }


    /**
     * @brief Allocate a Field to the advected function.
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_rp
     *      A DConstVectorFieldRTheta containing the values of the advection field
     *      on the logical index range axis.
     * @param [in] advection_field_xy_center
     *      A CoordXY containing the value of the advection field on the 
     *      physical index range axis at the O-point. 
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<R, Theta> advection_field_rp,
            CoordXY const& advection_field_xy_center,
            double dt) const
    {
        IdxRangeRTheta grid(get_idx_range<GridR, GridTheta>(allfdistribu));

        const int npoints_p = IdxRangeTheta(grid).size();
        IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
        IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_p)));


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemRTheta<X, Y> advection_field_xy(grid);

        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irp) {
            CoordRTheta const coord_rp(ddc::coordinate(irp));

            std::array<std::array<double, 2>, 2> J; // Jacobian matrix
            m_mapping.jacobian_matrix(coord_rp, J);
            std::array<std::array<double, 2>, 2> G; // Metric tensor
            m_mapping.metric_tensor(coord_rp, G);

            ddcHelper::get<X>(advection_field_xy)(irp)
                    = ddcHelper::get<R>(advection_field_rp)(irp) * J[1][1] / std::sqrt(G[1][1])
                      + ddcHelper::get<Theta>(advection_field_rp)(irp) * -J[1][0]
                                / std::sqrt(G[0][0]);
            ddcHelper::get<Y>(advection_field_xy)(irp)
                    = ddcHelper::get<R>(advection_field_rp)(irp) * -J[0][1] / std::sqrt(G[1][1])
                      + ddcHelper::get<Theta>(advection_field_rp)(irp) * J[0][0]
                                / std::sqrt(G[0][0]);
        });

        ddc::for_each(Opoint_grid, [&](IdxRTheta const irp) {
            ddcHelper::get<X>(advection_field_xy)(irp) = CoordX(advection_field_xy_center);
            ddcHelper::get<Y>(advection_field_xy)(irp) = CoordY(advection_field_xy_center);
        });

        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRTheta> const interpolator_ptr = m_interpolator.preallocate();

        // Initialise the feet
        FieldMemRTheta<CoordRTheta> feet_rp(grid);
        ddc::for_each(grid, [&](IdxRTheta const irp) { feet_rp(irp) = ddc::coordinate(irp); });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(get_field(feet_rp), advection_field_xy, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, get_const_field(feet_rp));

        return allfdistribu;
    }
};
