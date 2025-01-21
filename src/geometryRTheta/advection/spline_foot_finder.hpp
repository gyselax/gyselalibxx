// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "ifoot_finder.hpp"
#include "vector_mapper.hpp"

/**
 * @brief A class to find the foot of the characteristics on the @f$ (r,\theta) @f$
 * plane.
 *
 * The natural advection domain is the physical domain,
 * where the studied equation is given.
 * However, not all the mappings used are analytically invertible and inverting
 * the Jacobian matrix of the mapping could be costly and could introduce numerical
 * errors. That is why, we also introduce a pseudo-Cartesian domain.
 *
 * More details can be found in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * @tparam TimeStepper
 *      A child class of ITimeStepper providing a time integration method.
 * @tparam LogicalToPhysicalMapping
 *      A mapping from the logical domain to the physical domain.
 * @tparam LogicalToPseudoPhysicalMapping
 *      A mapping from the logical domain to the domain where the advection is
 *      carried out. This may be a pseudo-physical domain or the physical domain
 *      itself.
 *
 * @see BslAdvectionRTheta
 */
template <class TimeStepper, class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping>
class SplineFootFinder : public IFootFinder
{
    static_assert(is_mapping_v<LogicalToPhysicalMapping>);
    static_assert(is_mapping_v<LogicalToPseudoPhysicalMapping>);
    static_assert(is_analytical_mapping_v<LogicalToPseudoPhysicalMapping>);

private:
    using PseudoPhysicalToLogicalMapping = inverse_mapping_t<LogicalToPseudoPhysicalMapping>;

private:
    /**
     * @brief Tag the first dimension in the advection index range.
     */
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    /**
     * @brief Tag the second dimension in the advection index range.
     */
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_adv = Coord<X_adv, Y_adv>;

    using PseudoCartesianToCircular = CartesianToCircular<X_adv, Y_adv, R, Theta>;
    using PseudoPhysicalToPhysicalMapping
            = CombinedMapping<LogicalToPhysicalMapping, PseudoCartesianToCircular>;

    TimeStepper const& m_time_stepper;

    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    PseudoPhysicalToPhysicalMapping m_pseudo_physical_to_physical;

    SplineRThetaBuilder_host const& m_builder_advection_field;
    SplineRThetaEvaluatorConstBound_host const& m_evaluator_advection_field;


public:
    /**
     * @brief Instantiate a time integration method for the advection
     * operator.
     *
     * @param[in] time_stepper
     *      The time integration method used to solve the characteristic
     *      equation (ITimeStepper).
     * @param[in] logical_to_physical_mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] logical_to_pseudo_physical_mapping
     *      The mapping from the logical domain to the pseudo-physical domain.
     * @param[in] builder_advection_field
     *      The spline builder which computes the spline representation
     *      of the advection field.
     * @param[in] evaluator_advection_field
     *      The B-splines evaluator to evaluate the advection field.
     * @param[in] epsilon
     *      @f$ \varepsilon @f$ parameter used for the linearization of the
     *      advection field around the central point.
     *
     * @see ITimeStepper
     */
    SplineFootFinder(
            TimeStepper const& time_stepper,
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical_mapping,
            SplineRThetaBuilder_host const& builder_advection_field,
            SplineRThetaEvaluatorConstBound_host const& evaluator_advection_field,
            double epsilon = 1e-12)
        : m_time_stepper(time_stepper)
        , m_logical_to_pseudo_physical(logical_to_pseudo_physical_mapping)
        , m_pseudo_physical_to_logical(logical_to_pseudo_physical_mapping.get_inverse_mapping())
        , m_pseudo_physical_to_physical(
                  logical_to_physical_mapping,
                  PseudoCartesianToCircular(),
                  epsilon)
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
                m_pseudo_physical_to_physical);

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

        IdxRangeRTheta const idx_range_rp = get_idx_range<GridR, GridTheta>(feet);

        CoordXY_adv coord_center(m_logical_to_pseudo_physical(CoordRTheta(0, 0)));

        // The function describing how the value(s) are updated using the derivative.
        std::function<
                void(host_t<FieldRTheta<CoordRTheta>>,
                     host_t<DConstVectorFieldRTheta<X_adv, Y_adv>>,
                     double)>
                update_function = [&](host_t<FieldRTheta<CoordRTheta>> feet,
                                      host_t<DConstVectorFieldRTheta<X_adv, Y_adv>> advection_field,
                                      double dt) {
                    // Compute the characteristic feet at t^n:
                    ddc::for_each(idx_range_rp, [&](IdxRTheta const irp) {
                        CoordRTheta const coord_rp(feet(irp));
                        CoordXY_adv const coord_xy = m_logical_to_pseudo_physical(coord_rp);

                        CoordXY_adv const feet_xy = coord_xy - dt * advection_field(irp);

                        if (norm_inf(feet_xy - coord_center) < 1e-15) {
                            feet(irp) = CoordRTheta(0, 0);
                        } else {
                            feet(irp) = m_pseudo_physical_to_logical(feet_xy);
                            ddc::select<Theta>(feet(irp)) = ddcHelper::restrict_to_idx_range(
                                    ddc::select<Theta>(feet(irp)),
                                    IdxRangeTheta(idx_range_rp));
                        }
                    });

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
    void unify_value_at_center_pt(Field<T, IdxRangeRTheta, Kokkos::HostSpace> values) const
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
