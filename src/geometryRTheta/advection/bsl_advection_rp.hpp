// SPDX-License-Identifier: MIT

#pragma once

#include <functional>

#include <sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>

#include <directional_tag.hpp>
#include <geometry.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "advection_domain.hpp"
#include "foot_finder.hpp"
#include "i_interpolator_2d_rp.hpp"
#include "iadvectionrp.hpp"
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
template <class AdvectionDomain, class TimeStepper>
class BslAdvectionRP : public IAdvectionRP
{
private:
    using RDimX_adv = typename AdvectionDomain::RDimX_adv;
    using RDimY_adv = typename AdvectionDomain::RDimY_adv;

    AdvectionDomain const& m_advection_domain;

    PreallocatableSplineInterpolatorRP const& m_interpolator;

    IFootFinder<TimeStepper, AdvectionDomain> m_find_feet;

    FieldRP<CoordRP> m_feet_coords_rp;



public:
    /**
     * @brief Instantiate an advection operator.
     *
     * @param[in] grid
     *      The domain on which the advected function is defined.
     * @param [in] advection_domain
     *      The advection domain where the characteristic feet are computed.
     *      The AdvectionDomain object contains the Mapping information.
     * @param [in] function_interpolator
     *       The polar interpolator to interpolate the function once the
     *      characteristic computed.
     * @param[in] advection_builder
     *      The B-splines builder used to define the B-splines coefficients
     *      of the advection field.
     * @param[in] advection_evaluator
     *      The B-splines evaluator used to evaluate the advection field.
     * @param [in] time_stepper
     *       The time integration method used for the computation of the
     *      characteristic feet.
     */
    BslAdvectionRP(
            IDomainRP const& grid,
            AdvectionDomain const& advection_domain,
            PreallocatableSplineInterpolatorRP const& function_interpolator,
            SplineRPBuilder const& advection_builder,
            SplineRPEvaluator& advection_evaluator,
            TimeStepper& time_stepper)
        : m_advection_domain(advection_domain)
        , m_interpolator(function_interpolator)
        , m_find_feet(time_stepper, advection_domain, grid, advection_builder, advection_evaluator)
        , m_feet_coords_rp(grid)
    {
    }

    ~BslAdvectionRP() override = default;


    /**
     * @brief Allocate a ChunkSpan to the advected function.
     *
     * @param [in, out] allfdistribu
     *      A ChunkSpan containing the values of the function we want to advect.
     * @param [in] advection_field
     *      A VectorDViewRP containing the values of the advection field
     *      in the physical domain.
     * @param [in] dt
     *      A time step used.
     *
     * @return A ChunkSpan to allfdistribu advected on the time step given.
     */
    DSpanRP operator()(
            DSpanRP allfdistribu,
            VectorDViewRP<RDimX, RDimY> advection_field,
            double const dt)
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRP> const interpolator_ptr = m_interpolator.preallocate();

        // Initialisation the feet
        ddc::for_each(advection_field.domain(), [&](IndexRP const irp) {
            m_feet_coords_rp(irp) = ddc::coordinate(irp);
        });


        // Set the attributes of the m_find_feet object and compute the advection field in the
        // advection domain. ------------------------------------------------------------------------
        m_find_feet.set_vector_field(advection_field);

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet.update_feet(m_feet_coords_rp.span_view(), dt);

        m_find_feet.is_unified(m_feet_coords_rp.span_view());

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, m_feet_coords_rp.span_cview());

        return allfdistribu;
    }


    /**
     * @brief Get the values of current characteristic feet.
     *
     * This function must be called once advection computed,
     * so once the operator() called.
     *
     * @return A ChunkSpan to the computed values of the characteristic feet.
     */
    SpanRP<CoordRP> get_feet()
    {
        return m_feet_coords_rp.span_view();
    };
};
