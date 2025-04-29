// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "ipolar_foot_finder.hpp"
#include "mapping_tools.hpp"
#include "vector_index_tools.hpp"
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
 * @tparam SplineRThetaBuilderAdvection
 *      A 2D SplineBuilder to construct a spline on a polar domain.
 * @tparam SplineRThetaEvaluatorAdvection
 *      A 2D SplineEvaluator to evaluate a spline on a polar domain.
 *      A boundary condition must be provided in case the foot of the characteristic
 *      is found outside the domain.
 *
 * @see BslAdvectionPolar
 */
template <
        class TimeStepper,
        class LogicalToPhysicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplinePolarFootFinder
    : public IPolarFootFinder<
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
              ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
              typename TimeStepper::IdxRange,
              typename SplineRThetaBuilderAdvection::memory_space>
{
    static_assert(is_mapping_v<LogicalToPhysicalMapping>);
    static_assert(is_mapping_v<LogicalToPseudoPhysicalMapping>);
    static_assert(is_analytical_mapping_v<LogicalToPseudoPhysicalMapping>);
    static_assert(std::is_same_v<
                  typename SplineRThetaBuilderAdvection::memory_space,
                  typename SplineRThetaEvaluatorAdvection::memory_space>);
    static_assert(is_accessible_v<
                  typename SplineRThetaBuilderAdvection::exec_space,
                  LogicalToPhysicalMapping>);
    static_assert(is_accessible_v<
                  typename SplineRThetaBuilderAdvection::exec_space,
                  LogicalToPseudoPhysicalMapping>);
    static_assert(ddc::in_tags_v<
                  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
                  ddc::to_type_seq_t<typename TimeStepper::IdxRange>>);
    static_assert(ddc::in_tags_v<
                  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
                  ddc::to_type_seq_t<typename TimeStepper::IdxRange>>);

    using base_type = IPolarFootFinder<
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
            ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
            typename TimeStepper::IdxRange,
            typename SplineRThetaBuilderAdvection::memory_space>;

public:
    using typename base_type::GridR;
    using typename base_type::GridTheta;
    using typename base_type::IdxRangeOperator;
    using typename base_type::memory_space;
    using typename base_type::R;
    using typename base_type::Theta;
    using typename base_type::VectorIndexSetAdvectionDims;

private:
    using PseudoPhysicalToLogicalMapping = inverse_mapping_t<LogicalToPseudoPhysicalMapping>;

private:
    using ExecSpace = typename SplineRThetaBuilderAdvection::exec_space;
    /**
     * @brief Tag the first dimension in the advection domain.
     */
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    /**
     * @brief Tag the second dimension in the advection domain.
     */
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_adv = typename LogicalToPseudoPhysicalMapping::CoordResult;

    using PseudoCartesianBasis = ddc::to_type_seq_t<CoordXY_adv>;

    using CoordRTheta = Coord<R, Theta>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;

    using PseudoCartesianToCircular = CartesianToCircular<X_adv, Y_adv, R, Theta>;
    using PseudoPhysicalToPhysicalMapping
            = CombinedMapping<LogicalToPhysicalMapping, PseudoCartesianToCircular>;

    using BSplinesR = typename SplineRThetaBuilderAdvection::bsplines_type1;
    using BSplinesTheta = typename SplineRThetaBuilderAdvection::bsplines_type2;

    using IdxRangeSplineBatched
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeOperator>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    TimeStepper const& m_time_stepper;

    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    PseudoPhysicalToPhysicalMapping m_pseudo_physical_to_physical;

    SplineRThetaBuilderAdvection const& m_builder_advection_field;
    SplineRThetaEvaluatorAdvection const& m_evaluator_advection_field;

public:
    /**
     * @brief The type of a field of (r, theta) coordinates at every grid point, saved
     * on a compatible memory space.
     */
    using CFieldFeet = Field<CoordRTheta, IdxRangeOperator, memory_space>;

    /**
     * @brief The type of a constant field of (r, theta) coordinates at every grid point,
     * saved on a compatible memory space.
     */
    using CConstFieldFeet = ConstField<CoordRTheta, IdxRangeOperator, memory_space>;

    /**
     * @brief The type of a vector field defined on the pseudo-Cartesian basis at every
     * grid point, saved on a compatible memory space.
     */
    using DVectorFieldAdvection
            = DVectorField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;

    /**
     * @brief The type of a constant vector field defined on the pseudo-Cartesian basis
     * at every grid point, saved on a compatible memory space.
     */
    using DVectorConstFieldAdvection
            = DVectorConstField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;

    /**
     * @brief The type of 2 batched splines representing the x and y components of a vector
     * on the polar plane on a compatible memory space.
     */
    using VectorSplineCoeffsMem
            = DVectorFieldMem<IdxRangeSplineBatched, PseudoCartesianBasis, memory_space>;

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
     *      @f$ \varepsilon @f$ parameter used for the linearisation of the
     *      advection field around the central point.
     *
     * @see ITimeStepper
     */
    SplinePolarFootFinder(
            TimeStepper const& time_stepper,
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical_mapping,
            SplineRThetaBuilderAdvection const& builder_advection_field,
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
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
    void operator()(
            CFieldFeet feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt) const final
    {
        VectorSplineCoeffsMem advection_field_in_adv_domain_coefs(
                get_spline_idx_range(m_builder_advection_field));

        // Compute the advection field in the advection domain.
        auto advection_field_in_adv_domain = create_geometry_mirror_view(
                ExecSpace(),
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
        std::function<void(DVectorFieldAdvection, CConstFieldFeet)> dy
                = [&](DVectorFieldAdvection updated_advection_field, CConstFieldFeet feet) {
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

        CoordXY_adv coord_centre(m_logical_to_pseudo_physical(CoordRTheta(0, 0)));
        LogicalToPseudoPhysicalMapping logical_to_pseudo_physical_proxy
                = m_logical_to_pseudo_physical;
        PseudoPhysicalToLogicalMapping pseudo_physical_to_logical_proxy
                = m_pseudo_physical_to_logical;

        IdxRangeOperator idx_range = get_idx_range(feet);
        IdxRangeTheta idx_range_theta(idx_range);

        // The function describing how the value(s) are updated using the derivative.
        std::function<void(CFieldFeet, DVectorConstFieldAdvection, double)> update_function
                = [&](CFieldFeet feet, DVectorConstFieldAdvection advection_field, double dt) {
                      // Compute the characteristic feet at t^n:
                      ddc::parallel_for_each(
                              ExecSpace(),
                              get_idx_range(feet),
                              KOKKOS_LAMBDA(IdxOperator const idx) {
                                  CoordRTheta const coord_rtheta(feet(idx));
                                  CoordXY_adv const coord_xy
                                          = logical_to_pseudo_physical_proxy(coord_rtheta);

                                  CoordXY_adv const feet_xy = coord_xy - dt * advection_field(idx);

                                  if (norm_inf(feet_xy - coord_centre) < 1e-15) {
                                      feet(idx) = CoordRTheta(0, 0);
                                  } else {
                                      feet(idx) = pseudo_physical_to_logical_proxy(feet_xy);
                                      ddc::select<Theta>(feet(idx))
                                              = ddcHelper::restrict_to_idx_range(
                                                      ddc::select<Theta>(feet(idx)),
                                                      idx_range_theta);
                                  }
                              });

                      // Treatment to conserve the C0 property of the advected function:
                      unify_value_at_centre_pt(feet);
                      // Test if the values are the same at the centre point
                      is_unified(feet);
                  };


        // Solve the characteristic equation
        m_time_stepper.update(ExecSpace(), feet, dt, dy, update_function);

        is_unified(feet);
    }



    /**
     * @brief Check if the values at the centre point are the same.
     *
     *  For polar geometry, to ensure continuity at the centre point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  This function check if for @f$ r= 0 @f$, the values @f$ \forall \theta @f$ are the same.
     *
     *  @param[in] values
     *      A table of values we want to check if the centre point has
     *      an unique value.
     *
     */
    template <class T>
    void is_unified(Field<T, IdxRangeRTheta, memory_space> const& values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        IdxR r0_idx = r_idx_range.front();
        if (Kokkos::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            ddc::parallel_for_each(
                    ExecSpace(),
                    theta_idx_range,
                    KOKKOS_LAMBDA(const IdxTheta itheta) {
                        if (norm_inf(
                                    values(r0_idx, itheta)
                                    - values(r0_idx, theta_idx_range.front()))
                            > 1e-15) {
                            Kokkos::printf("WARNING ! -> Discontinuous at the centre point.");
                        }
                        KOKKOS_ASSERT(
                                values(r0_idx, itheta) == values(r0_idx, theta_idx_range.front()));
                    });
        }
    }


    /**
     * @brief Replace the value at @f$  (r=0, \theta)@f$  point
     *  by the value at @f$ (r=0,0) @f$ for all @f$ \theta @f$.
     *
     *  For polar geometry, to ensure continuity at the centre point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  As the computation of the values of a table can induces machine errors,
     *  this function is useful to reset the values at the central point at
     *  the same value.
     *
     *  @param[in, out] values
     *      The table of values we want to unify at the central point.
     */
    template <class T>
    void unify_value_at_centre_pt(Field<T, IdxRangeRTheta, memory_space> values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        IdxR r0_idx = r_idx_range.front();
        if (std::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            ddc::parallel_for_each(
                    ExecSpace(),
                    theta_idx_range,
                    KOKKOS_LAMBDA(const IdxTheta itheta) {
                        values(r0_idx, itheta) = values(r0_idx, theta_idx_range.front());
                    });
        }
    }
};
