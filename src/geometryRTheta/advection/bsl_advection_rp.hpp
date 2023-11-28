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
template <class FootFinder>
class BslAdvectionRP : public IAdvectionRP
{
private:
    PreallocatableSplineInterpolatorRP const& m_interpolator;

    FootFinder const& m_find_feet;


public:
    /**
     * @brief Instantiate an advection operator.
     *
     * @param [in] function_interpolator
     *       The polar interpolator to interpolate the function once the
     *      characteristic computed.
     * @param[in] foot_finder
     *      An IFootFinder which computes the characteristic feet.
     *
     * @tparam IFootFinder
     *      A child class of IFootFinder.
     */
    BslAdvectionRP(
            PreallocatableSplineInterpolatorRP const& function_interpolator,
            FootFinder const& foot_finder)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
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
    DSpanRP operator()(DSpanRP allfdistribu, VectorDViewRP<RDimX, RDimY> advection_field, double dt)
            const
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolatorRP> const interpolator_ptr = m_interpolator.preallocate();

        // Initialisation the feet
        FieldRP<CoordRP> feet_rp(advection_field.domain());
        ddc::for_each(advection_field.domain(), [&](IndexRP const irp) {
            feet_rp(irp) = ddc::coordinate(irp);
        });

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(feet_rp.span_view(), advection_field, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        (*interpolator_ptr)(allfdistribu, feet_rp.span_cview());

        return allfdistribu;
    }
};
