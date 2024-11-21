// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include <sll/mapping/cartesian_to_pseudo_cartesian.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "ifoot_finder.hpp"
#include "vector_mapper.hpp"

/**
 * @brief Define a base class for all the time integration methods used for the advection.
 *
 * @tparam TimeStepper
 *      A child class of ITimeStepper providing a time integration method.
 * @tparam AdvectionDomain
 *      A child class of AdvectionDomain providing the informations about the advection index range.
 *
 * @see BslAdvectionRTheta
 */
template <class TimeStepper, class AdvectionDomain, class Mapping>
class SplineFootFinder : public IFootFinder
{
private:
    /**
     * @brief Tag the first dimension in the advection index range.
     */
    using X_adv = typename AdvectionDomain::X_adv;
    /**
     * @brief Tag the second dimension in the advection index range.
     */
    using Y_adv = typename AdvectionDomain::Y_adv;

    using CircularToPseudoCartesian = CircularToCartesian<X_adv, Y_adv, R, Theta>;
    using AdvectionMapping = CartesianToPseudoCartesian<Mapping, CircularToPseudoCartesian>;

    TimeStepper const& m_time_stepper;

    AdvectionDomain const& m_advection_domain;
    AdvectionMapping m_mapping;

    SplineRThetaBuilder const& m_builder_advection_field;
    SplineRThetaEvaluatorConstBound const& m_evaluator_advection_field;


public:
    /**
     * @brief Instantiate a time integration method for the advection
     * operator.
     *
     * @param[in] time_stepper
     *      The time integration method used to solve the characteristic
     *      equation (ITimeStepper).
     * @param[in] advection_domain
     *      An AdvectionDomain object which defines in which index range we
     *      advect the characteristics.
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] builder_advection_field
     *      The spline builder which computes the spline representation
     *      of the advection field.
     * @param[in] evaluator_advection_field
     *      The B-splines evaluator to evaluate the advection field.
     * @param[in] epsilon
     *      @f$ \varepsilon @f$ parameter used for the linearization of the
     *      advection field around the central point.
     *
     * @tparam TimeStepper
     *      A child class of ITimeStepper providing a time integration method.
     * @tparam AdvectionDomain
     *      A child class of AdvectionDomain providing the informations about the advection index range.
     *
     * @see ITimeStepper
     */
    SplineFootFinder(
            TimeStepper const& time_stepper,
            AdvectionDomain const& advection_domain,
            Mapping const& mapping,
            SplineRThetaBuilder const& builder_advection_field,
            SplineRThetaEvaluatorConstBound const& evaluator_advection_field,
            double epsilon = 1e-12)
        : m_time_stepper(time_stepper)
        , m_advection_domain(advection_domain)
        , m_mapping(CircularToPseudoCartesian(), mapping, epsilon)
        , m_builder_advection_field(builder_advection_field)
        , m_evaluator_advection_field(evaluator_advection_field)
    {
    }

    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * From the advection field in the physical index range, compute the advection field
     * in the right index range an compute its B-splines coefficients.
     * Then, use the given time integration method (time_stepper) to solve the
     * characteristic equation over @f$ dt @f$.
     *
     * @param[in, out] feet
     *      On input: the mesh points.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the physical index range.
     * @param[in] dt
     *      The time step.
     */
    void operator()(
            host_t<FieldRTheta<CoordRTheta>> feet,
            host_t<DConstVectorFieldRTheta<X, Y>> advection_field,
            double dt) const final
    {
        host_t<VectorSplineCoeffsMem2D<X_adv, Y_adv>> advection_field_in_adv_domain_coefs(
                get_spline_idx_range(m_builder_advection_field));

        // Compute the advection field in the advection domain.
        auto advection_field_in_adv_domain = create_geometry_mirror_view(
                Kokkos::DefaultHostExecutionSpace(),
                advection_field,
                m_mapping);

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<X_adv>(advection_field_in_adv_domain_coefs),
                ddcHelper::get<X_adv>(get_const_field(advection_field_in_adv_domain)));
        m_builder_advection_field(
                ddcHelper::get<Y_adv>(advection_field_in_adv_domain_coefs),
                ddcHelper::get<Y_adv>(get_const_field(advection_field_in_adv_domain)));


        // The function describing how the derivative of the evolve function is calculated.
        std::function<
                void(host_t<DVectorFieldRTheta<X_adv, Y_adv>>,
                     host_t<ConstFieldRTheta<CoordRTheta>>)>
                dy = [&](host_t<DVectorFieldRTheta<X_adv, Y_adv>> updated_advection_field,
                         host_t<ConstFieldRTheta<CoordRTheta>> feet) {
                    m_evaluator_advection_field(
                            get_field(ddcHelper::get<X_adv>(updated_advection_field)),
                            get_const_field(feet),
                            get_const_field(
                                    ddcHelper::get<X_adv>(advection_field_in_adv_domain_coefs)));
                    m_evaluator_advection_field(
                            get_field(ddcHelper::get<Y_adv>(updated_advection_field)),
                            get_const_field(feet),
                            get_const_field(
                                    ddcHelper::get<Y_adv>(advection_field_in_adv_domain_coefs)));
                };

        // The function describing how the value(s) are updated using the derivative.
        std::function<
                void(host_t<FieldRTheta<CoordRTheta>>,
                     host_t<DConstVectorFieldRTheta<X_adv, Y_adv>>,
                     double)>
                update_function = [&](host_t<FieldRTheta<CoordRTheta>> feet,
                                      host_t<DConstVectorFieldRTheta<X_adv, Y_adv>> advection_field,
                                      double dt) {
                    // Compute the characteristic feet at t^n:
                    m_advection_domain.advect_feet(feet, advection_field, dt);

                    // Treatment to conserve the C0 property of the advected function:
                    unify_value_at_center_pt(feet);
                    // Test if the values are the same at the center point
                    is_unified(feet);
                };


        // Solve the characteristic equation
        m_time_stepper.update(Kokkos::DefaultHostExecutionSpace(), feet, dt, dy, update_function);

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
    void is_unified(Field<T, IdxRangeRTheta, Kokkos::HostSpace> const& values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
            ddc::for_each(theta_idx_range, [&](const IdxTheta ip) {
                if (norm_inf(
                            values(r_idx_range.front(), ip)
                            - values(r_idx_range.front(), theta_idx_range.front()))
                    > 1e-15) {
                    std::cout << "WARNING ! -> Discontinous at the center point." << std::endl;
                }
                assert(values(r_idx_range.front(), ip)
                       == values(r_idx_range.front(), theta_idx_range.front()));
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
    void unify_value_at_center_pt(FieldRTheta<T> values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
            ddc::for_each(theta_idx_range, [&](const IdxTheta ip) {
                values(r_idx_range.front(), ip)
                        = values(r_idx_range.front(), theta_idx_range.front());
            });
        }
    }
};
