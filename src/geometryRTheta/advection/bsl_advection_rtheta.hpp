// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "i_interpolator_rtheta.hpp"
#include "iadvection_rtheta.hpp"
#include "metric_tensor_evaluator.hpp"
#include "spline_interpolator_rtheta.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



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
 * The feet can be advected on different domains.
 * Theses domains are defined in the AdvectionDomain class.
 *
 * The interpolation of the function is always done in the logical index range.
 *
 *
 *
 * @see IFootFinder
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
     *      The mapping function from the logical domain to the physical
     *      domain. 
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
     *      on the physical domain axes.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<X, Y> advection_field_xy,
            double dt) const override
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRTheta> const interpolator_ptr = m_interpolator.preallocate();

        // Initialise the feet
        FieldMemRTheta<CoordRTheta> feet_rtheta_alloc(get_idx_range(advection_field_xy));
        FieldRTheta<CoordRTheta> feet_rtheta = get_field(feet_rtheta_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(advection_field_xy),
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    feet_rtheta(irtheta) = ddc::coordinate(irtheta);
                });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(feet_rtheta, get_const_field(advection_field_xy), dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(get_field(allfdistribu), get_const_field(feet_rtheta));

        return allfdistribu;
    }


    /**
     * @brief Allocate a Field to the advected function.
     *
     * @param [in, out] allfdistribu_host
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_rtheta
     *      A DConstVectorFieldRTheta containing the values of the advection field
     *      on the logical index range axis.
     * @param [in] advection_field_xy_centre
     *      A CoordXY containing the value of the advection field on the 
     *      physical index range axis at the O-point. 
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu_host,
            host_t<DConstVectorFieldRTheta<R, Theta>> advection_field_rtheta,
            CoordXY const& advection_field_xy_centre,
            double dt) const override
    {
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeRTheta grid(get_idx_range<GridR, GridTheta>(allfdistribu_host));

        const int npoints_theta = IdxRangeTheta(grid).size();
        IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
        IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_theta)));


        // Convert advection field on RTheta to advection field on XY
        host_t<DVectorFieldMemRTheta<X, Y>> advection_field_xy_host(grid);

        InverseJacobianMatrix<Mapping, CoordRTheta> inv_jacobian_matrix(m_mapping);

        ddc::for_each(grid_without_Opoint, [&](IdxRTheta const irtheta) {
            CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

            std::array<std::array<double, 2>, 2> inv_J = inv_jacobian_matrix(coord_rtheta);
            double const jacobian = m_mapping.jacobian(coord_rtheta);

            ddcHelper::get<X>(advection_field_xy_host)(irtheta)
                    = ddcHelper::get<R>(advection_field_rtheta)(irtheta) * inv_J[0][0] * jacobian
                      + ddcHelper::get<Theta>(advection_field_rtheta)(irtheta) * inv_J[1][0]
                                * jacobian;
            ddcHelper::get<Y>(advection_field_xy_host)(irtheta)
                    = ddcHelper::get<R>(advection_field_rtheta)(irtheta) * inv_J[0][1] * jacobian
                      + ddcHelper::get<Theta>(advection_field_rtheta)(irtheta) * inv_J[1][1]
                                * jacobian;
        });

        ddc::for_each(Opoint_grid, [&](IdxRTheta const irtheta) {
            ddcHelper::get<X>(advection_field_xy_host)(irtheta) = CoordX(advection_field_xy_centre);
            ddcHelper::get<Y>(advection_field_xy_host)(irtheta) = CoordY(advection_field_xy_centre);
        });

        auto allfdistribu = ddc::
                create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), allfdistribu_host);

        auto advection_field_xy = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(advection_field_xy_host));

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        ddc::parallel_deepcopy(allfdistribu_host, get_const_field(allfdistribu));

        Kokkos::Profiling::popRegion();

        return allfdistribu_host;
    }
};
