#pragma once

#include <functional>

#include "ifoot_finder.hpp"

/**
 * @brief Define a base class for all the time integration methods used for the advection.
 *
 * @tparam TimeStepper
 *      A child class of ITimeStepper providing a time integration method.
 * @tparam AdvectionDomain
 *      A child class of AdvectionDomain providing the informations about the advection domain.
 *
 * @see BslAdvectionRP
 */
template <class TimeStepper, class AdvectionDomain>
class SplineFootFinder : public IFootFinder
{
private:
    /**
     * @brief Tag the first dimension in the advection domain.
     */
    using RDimX_adv = typename AdvectionDomain::RDimX_adv;
    /**
     * @brief Tag the second dimension in the advection domain.
     */
    using RDimY_adv = typename AdvectionDomain::RDimY_adv;


    TimeStepper const& m_time_stepper;

    AdvectionDomain const& m_advection_domain;

    SplineRPBuilder const& m_builder_advection_field;
    SplineRPEvaluator const& m_evaluator_advection_field;


public:
    /**
     * @brief Instantiate a time integration method for the advection
     * operator.
     *
     * @param[in] time_stepper
     *      The time integration method used to solve the characteristic
     *      equation (ITimeStepper).
     * @param[in] advection_domain
     *      An AdvectionDomain object which defines in which domain we
     *      advect the characteristics.
     * @param[in] builder_advection_field
     *      The spline builder which computes the spline representation
     *      of the advection field.
     * @param[in] evaluator_advection_field
     *      The B-splines evaluator to evaluate the advection field.
     *
     * @tparam TimeStepper
     *      A child class of ITimeStepper providing a time integration method.
     * @tparam AdvectionDomain
     *      A child class of AdvectionDomain providing the informations about the advection domain.
     *
     * @see ITimeStepper
     */
    SplineFootFinder(
            TimeStepper const& time_stepper,
            AdvectionDomain const& advection_domain,
            SplineRPBuilder const& builder_advection_field,
            SplineRPEvaluator const& evaluator_advection_field)
        : m_time_stepper(time_stepper)
        , m_advection_domain(advection_domain)
        , m_builder_advection_field(builder_advection_field)
        , m_evaluator_advection_field(evaluator_advection_field)
    {
    }

    ~SplineFootFinder() {};


    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * From the advection field in the physical domain, compute the advection field
     * in the right domain an compute its B-splines coefficients.
     * Then, use the given time integration method (time_stepper) to solve the
     * characteristic equation over @f$ dt @f$.
     *
     * @param[in, out] feet
     *      On input: the mesh points.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the physical domain.
     * @param[in] dt
     *      The time step.
     */
    void operator()(SpanRP<CoordRP> feet, VectorDViewRP<RDimX, RDimY> advection_field, double dt)
            const final
    {
        VectorDFieldRP<RDimX_adv, RDimY_adv> advection_field_in_adv_dom(advection_field.domain());
        VectorSpline2D<RDimX_adv, RDimY_adv> advection_field_in_adv_dom_coefs(
                m_builder_advection_field.spline_domain());

        // Compute the advection field in the advection domain.
        m_advection_domain
                .compute_advection_field(advection_field, advection_field_in_adv_dom.span_view());

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<RDimX_adv>(advection_field_in_adv_dom_coefs),
                ddcHelper::get<RDimX_adv>(advection_field_in_adv_dom));
        m_builder_advection_field(
                ddcHelper::get<RDimY_adv>(advection_field_in_adv_dom_coefs),
                ddcHelper::get<RDimY_adv>(advection_field_in_adv_dom));


        // The function describing how the derivative of the evolve function is calculated.
        std::function<void(VectorDSpanRP<RDimX_adv, RDimY_adv>, ViewRP<CoordRP>)> dy
                = [&](VectorDSpanRP<RDimX_adv, RDimY_adv> updated_advection_field,
                      ViewRP<CoordRP> feet) {
                      m_evaluator_advection_field(
                              ddcHelper::get<RDimX_adv>(updated_advection_field).span_view(),
                              feet.span_cview(),
                              ddcHelper::get<RDimX_adv>(advection_field_in_adv_dom_coefs)
                                      .span_cview());
                      m_evaluator_advection_field(
                              ddcHelper::get<RDimY_adv>(updated_advection_field).span_view(),
                              feet.span_cview(),
                              ddcHelper::get<RDimY_adv>(advection_field_in_adv_dom_coefs)
                                      .span_cview());
                  };

        // The function describing how the value(s) are updated using the derivative.
        std::function<void(SpanRP<CoordRP>, VectorDViewRP<RDimX_adv, RDimY_adv>, double)>
                update_function = [&](SpanRP<CoordRP> feet,
                                      VectorDViewRP<RDimX_adv, RDimY_adv> advection_field,
                                      double dt) {
                    // Compute the characteristic feet at t^n:
                    m_advection_domain.advect_feet(feet, advection_field, dt);

                    // Treatment to conserve the C0 property of the advected function:
                    unify_value_at_center_pt(feet);
                    // Test if the values are the same at the center point
                    is_unified(feet);
                };


        // Solve the characteristic equation
        m_time_stepper.update(Kokkos::Serial(), feet, dt, dy, update_function);

        is_unified(feet);
    }



private:
    /**
     * @brief Check if the values at the center point are the same.
     *
     *  For polar geometry, to ensure continuity at the center point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  This function check if for @f$ r= 0 @f$, the values @f$ \forall \theta @f$ are the same.
     *
     *  @param[in] values
     *      A table of values we want to check if the center point has
     *      an unique value.
     *
     */
    template <class T>
    void is_unified(SpanRP<T> const& values) const
    {
        auto const r_domain = ddc::get_domain<IDimR>(values);
        auto const theta_domain = ddc::get_domain<IDimP>(values);
        if (std::fabs(ddc::coordinate(r_domain.front())) < 1e-15) {
            ddc::for_each(theta_domain, [&](const IndexP ip) {
                if (norm_inf(
                            values(r_domain.front(), ip)
                            - values(r_domain.front(), theta_domain.front()))
                    > 1e-15) {
                    std::cout << "WARNING ! -> Discontinous at the center point." << std::endl;
                }
                assert(values(r_domain.front(), ip)
                       == values(r_domain.front(), theta_domain.front()));
            });
        }
    }


    /**
     * @brief Replace the value at @f$  (r=0, \theta)@f$  point
     *  by the value at @f$ (r=0,0) @f$ for all @f$ \theta @f$.
     *
     *  For polar geometry, to ensure continuity at the center point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  As the computation of the values of a table can induces machine errors,
     *  this function is useful to reset the values at the central point at
     *  the same value.
     *
     *  @param[in, out] values
     *      The table of values we want to unify at the central point.
     */
    template <class T>
    void unify_value_at_center_pt(SpanRP<T> values) const
    {
        auto const r_domain = ddc::get_domain<IDimR>(values);
        auto const theta_domain = ddc::get_domain<IDimP>(values);
        if (std::fabs(ddc::coordinate(r_domain.front())) < 1e-15) {
            ddc::for_each(theta_domain, [&](const IndexP ip) {
                values(r_domain.front(), ip) = values(r_domain.front(), theta_domain.front());
            });
        }
    }
};
