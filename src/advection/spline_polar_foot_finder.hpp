// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "ipolar_foot_finder.hpp"
#include "itimestepper.hpp"
#include "l_norm_tools.hpp"
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
 * @tparam TimeStepperBuilder
 *      A time stepper builder indicating which time integration method should be
 *      applied to solve the characteristic equation. 
 * @tparam PseudoPhysicalToAdvectionDomainMapping
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
        class IdxRangeBatched,
        class TimeStepperBuilder,
        concepts::Mapping PseudoPhysicalToAdvectionDomainMapping,
        concepts::AnalyticalMapping LogicalToPseudoPhysicalMapping,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplinePolarFootFinder
    : public IPolarFootFinder<
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
              ddc::to_type_seq_t<typename PseudoPhysicalToAdvectionDomainMapping::CoordResult>,
              IdxRangeBatched,
              typename SplineRThetaBuilderAdvection::memory_space>
{
    static_assert(is_timestepper_builder_v<TimeStepperBuilder>);
    static_assert(std::is_same_v<
                  typename SplineRThetaBuilderAdvection::memory_space,
                  typename SplineRThetaEvaluatorAdvection::memory_space>);
    static_assert(is_accessible_v<
                  typename SplineRThetaBuilderAdvection::exec_space,
                  PseudoPhysicalToAdvectionDomainMapping>);
    static_assert(is_accessible_v<
                  typename SplineRThetaBuilderAdvection::exec_space,
                  LogicalToPseudoPhysicalMapping>);
    static_assert(ddc::in_tags_v<
                  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
                  ddc::to_type_seq_t<IdxRangeBatched>>);
    static_assert(ddc::in_tags_v<
                  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
                  ddc::to_type_seq_t<IdxRangeBatched>>);
    static_assert(
            SplineRThetaBuilderAdvection::builder_type1::s_nbe_xmin == 0,
            "This class is designed to work with a spline builder which does not require "
            "additional information at the boundaries (e.g. Hermite boundary conditions require "
            "information about the derivatives and therefore will not work with this class. Please "
            "check the choice of boundary conditions).");
    static_assert(
            SplineRThetaBuilderAdvection::builder_type1::s_nbe_xmax == 0,
            "This class is designed to work with a spline builder which does not require "
            "additional information at the boundaries (e.g. Hermite boundary conditions require "
            "information about the derivatives and therefore will not work with this class. Please "
            "check the choice of boundary conditions).");
    static_assert(
            SplineRThetaBuilderAdvection::builder_type1::s_bc_xmin != ddc::BoundCond::PERIODIC,
            "Periodic boundary conditions in the radial direction are nonsensical.");
    static_assert(
            SplineRThetaBuilderAdvection::builder_type1::s_bc_xmax != ddc::BoundCond::PERIODIC,
            "Periodic boundary conditions in the radial direction are nonsensical.");
    static_assert(
            SplineRThetaBuilderAdvection::builder_type2::s_bc_xmin == ddc::BoundCond::PERIODIC,
            "Expected periodic boundary conditions in the poloidal direction.");
    static_assert(
            SplineRThetaBuilderAdvection::builder_type2::s_bc_xmax == ddc::BoundCond::PERIODIC,
            "Expected periodic boundary conditions in the poloidal direction.");

    using base_type = IPolarFootFinder<
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
            ddc::to_type_seq_t<typename PseudoPhysicalToAdvectionDomainMapping::CoordResult>,
            IdxRangeBatched,
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

public:
    /// @brief Execution space.
    using ExecSpace = typename SplineRThetaBuilderAdvection::exec_space;

private:
    using MemSpace = typename ExecSpace::memory_space;
    /**
     * @brief Tag the first dimension in the advection domain.
     */
    using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
    /**
     * @brief Tag the second dimension in the advection domain.
     */
    using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;
    /**
     * @brief The coordinate type associated to the dimensions in the advection domain.
     */
    using CoordXY_pc = typename LogicalToPseudoPhysicalMapping::CoordResult;

    using CoordRTheta = Coord<R, Theta>;

    using PolarBasis = ddc::to_type_seq_t<CoordRTheta>;
    using PseudoCartesianBasis = ddc::to_type_seq_t<CoordXY_pc>;

    // Local X_pc/Y_pc is equal to global X_pc, Y_pc for non-invertible mappings or
    // directly equal to X,Y the physical coordinates
    using X_pc = ddc::type_seq_element_t<0, PseudoCartesianBasis>;
    using Y_pc = ddc::type_seq_element_t<1, PseudoCartesianBasis>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;

    using PseudoCartesianToCircular = CartesianToCircular<X_pc, Y_pc, R, Theta>;

    using BSplinesR = typename SplineRThetaBuilderAdvection::bsplines_type1;
    using BSplinesTheta = typename SplineRThetaBuilderAdvection::bsplines_type2;

    using IdxRangeSplineBatched
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeOperator>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    using TimeStepper = typename TimeStepperBuilder::
            template time_stepper_t<CoordRTheta, DVector<AdvDim1, AdvDim2>>;

    TimeStepperBuilder const& m_time_stepper_builder;

    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    PseudoPhysicalToAdvectionDomainMapping m_pseudo_physical_to_adv_domain;

    SplineRThetaBuilderAdvection const& m_builder_advection_field;
    SplineRThetaEvaluatorAdvection const& m_evaluator_advection_field;

public:
    /**
     * @brief The type of a field of (r, theta) coordinates at every grid point, saved
     * on a compatible memory space.
     */
    using CFieldFeet = Field<CoordRTheta, IdxRangeOperator, memory_space>;

    using VectorFieldBatchedSplineCoefMem
            = DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>;
    using VectorFieldBatchedSplineCoef
            = DVectorField<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>;
    using VectorConstFieldSplineCoef = DVectorConstField<
            IdxRange<BSplinesR, BSplinesTheta>,
            VectorIndexSetAdvectionDims,
            memory_space>;

public:
    /**
     * @brief Instantiate a time integration method for the advection
     * operator.
     *
     * @param[in] idx_range_operator
     *      The index range on which the operator should act.
     * @param[in] time_stepper_builder
     *      A builder for the time integration method used for the
     *      characteristic equation. 
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
            IdxRangeBatched const& idx_range_operator,
            TimeStepperBuilder const& time_stepper_builder,
            PseudoPhysicalToAdvectionDomainMapping const&
                    pseudo_physical_to_advection_domain_mapping,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical_mapping,
            SplineRThetaBuilderAdvection const& builder_advection_field,
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field)
        : m_time_stepper_builder(time_stepper_builder)
        , m_logical_to_pseudo_physical(logical_to_pseudo_physical_mapping)
        , m_pseudo_physical_to_logical(logical_to_pseudo_physical_mapping.get_inverse_mapping())
        , m_pseudo_physical_to_adv_domain(pseudo_physical_to_advection_domain_mapping)
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
        static_assert(ddc::type_seq_size_v<VectorIndexSetAdvectionDims> == 2);
        using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
        using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

        VectorFieldBatchedSplineCoefMem advection_field_coefs_alloc(
                m_builder_advection_field.batched_spline_domain(get_idx_range(advection_field)));
        VectorFieldBatchedSplineCoef advection_field_coefs(advection_field_coefs_alloc);

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));

        CoordXY_pc coord_centre(m_logical_to_pseudo_physical(CoordRTheta(0, 0)));
        LogicalToPseudoPhysicalMapping logical_to_pseudo_physical_proxy
                = m_logical_to_pseudo_physical;
        PseudoPhysicalToLogicalMapping pseudo_physical_to_logical_proxy
                = m_pseudo_physical_to_logical;

        IdxRangeOperator idx_range = get_idx_range(feet);
        IdxRangeTheta idx_range_theta(idx_range);

        TimeStepper time_stepper = m_time_stepper_builder.template preallocate<TimeStepper>();

        SplineRThetaEvaluatorAdvection const& evaluator_advection_field_proxy
                = m_evaluator_advection_field;
        PseudoPhysicalToAdvectionDomainMapping pseudo_physical_to_adv_domain
                = m_pseudo_physical_to_adv_domain;

        // Compute the characteristic feet at t^n:
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) {
                    IdxBatch idx_batch(idx);
                    IdxRTheta idx_rtheta(idx);
                    VectorConstFieldSplineCoef advection_field_coefs_slice
                            = get_const_field(advection_field_coefs)[idx_batch];
                    // The function describing how the derivative of the evolve function is calculated.
                    auto dy = [&](DVector<AdvDim1, AdvDim2>& updated_advection_field,
                                  CoordRTheta const& foot) {
                        ddcHelper::get<AdvDim1>(updated_advection_field)
                                = evaluator_advection_field_proxy(
                                        foot,
                                        ddcHelper::get<AdvDim1>(advection_field_coefs_slice));
                        ddcHelper::get<AdvDim2>(updated_advection_field)
                                = evaluator_advection_field_proxy(
                                        foot,
                                        ddcHelper::get<AdvDim2>(advection_field_coefs_slice));
                    };

                    // The function describing how the value(s) are updated using the derivative.
                    auto update_function = [&](CoordRTheta& foot_rtheta,
                                               DVector<AdvDim1, AdvDim2> const& advection_field,
                                               double dt) {
                        double radial_coord = ddc::select<R>(foot_rtheta);
                        CoordXY_pc foot_xy;
                        constexpr bool adv_provided_on_polar
                                = (ddc::in_tags_v<AdvDim1, PolarBasis>)&&(
                                        ddc::in_tags_v<AdvDim2, PolarBasis>);
                        if (!adv_provided_on_polar or radial_coord > 1e-15) {
                            // Ensure coord is inside the domain as splines can't extrapolate
                            // derivates (clamping)
                            CoordRTheta advection_location_for_mapping(
                                    Kokkos::
                                            min(ddc::select<R>(foot_rtheta),
                                                ddc::discrete_space<BSplinesR>().rmax()),
                                    ddc::select<Theta>(foot_rtheta));
                            using CoordJ =
                                    typename PseudoPhysicalToAdvectionDomainMapping::CoordJacobian;
                            DVector<X_pc, Y_pc> advection_field_xy;
                            if constexpr (std::is_same_v<CoordRTheta, CoordJ>) {
                                advection_field_xy = to_vector_space<VectorIndexSet<X_pc, Y_pc>>(
                                        pseudo_physical_to_adv_domain,
                                        advection_location_for_mapping,
                                        advection_field);
                            }
                            CoordXY_pc const coord_xy
                                    = logical_to_pseudo_physical_proxy(foot_rtheta);
                            if constexpr (!std::is_same_v<CoordRTheta, CoordJ>) {
                                advection_field_xy = to_vector_space<VectorIndexSet<X_pc, Y_pc>>(
                                        pseudo_physical_to_adv_domain,
                                        coord_xy,
                                        advection_field);
                            }
                            foot_xy = coord_xy - dt * advection_field_xy;
                        } else {
                            // TODO
                        }

                        if (norm_inf(foot_xy - coord_centre) < 1e-15) {
                            foot_rtheta = CoordRTheta(0, 0);
                        } else {
                            foot_rtheta = pseudo_physical_to_logical_proxy(foot_xy);
                            ddc::select<Theta>(foot_rtheta) = ddcHelper::restrict_to_idx_range(
                                    ddc::select<Theta>(foot_rtheta),
                                    idx_range_theta);
                        }
                    };

                    feet(idx) = ddc::coordinate(idx_rtheta);
                    // Solve the characteristic equation
                    time_stepper.update(feet(idx), dt, dy, update_function);
                });

        // Treatment to conserve the C0 property of the advected function:
        unify_value_at_centre_pt(feet);
        // Test if the values are the same at the centre point
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
    static void is_unified(Field<T, IdxRangeOperator, memory_space> const& values)
    {
        IdxRangeOperator full_idx_range = get_idx_range(values);
        IdxRangeBatch const batched_idx_range(full_idx_range);
        IdxRangeR const r_idx_range(full_idx_range);
        IdxRangeTheta const theta_idx_range(full_idx_range);
        IdxR r0_idx = r_idx_range.front();
        IdxTheta theta0_idx = theta_idx_range.front();
        if (Kokkos::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    ExecSpace(),
                    batched_idx_range,
                    KOKKOS_LAMBDA(const IdxBatch ib) {
                        for (IdxTheta itheta : theta_idx_range) {
                            if (norm_inf(
                                        values(ib, r0_idx, itheta) - values(ib, r0_idx, theta0_idx))
                                > 1e-15) {
                                Kokkos::printf("WARNING ! -> Discontinuous at the centre point.");
                            }
                            KOKKOS_ASSERT(
                                    values(ib, r0_idx, itheta) == values(ib, r0_idx, theta0_idx));
                        }
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
    static void unify_value_at_centre_pt(Field<T, IdxRangeOperator, memory_space> values)
    {
        IdxRangeOperator full_idx_range = get_idx_range(values);
        IdxRangeBatch const batched_idx_range(full_idx_range);
        IdxRangeR const r_idx_range(full_idx_range);
        IdxRangeTheta const theta_idx_range(full_idx_range);
        IdxR r0_idx = r_idx_range.front();
        IdxTheta theta0_idx = theta_idx_range.front();
        if (std::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    ExecSpace(),
                    batched_idx_range,
                    KOKKOS_LAMBDA(const IdxBatch ib) {
                        for (IdxTheta itheta : theta_idx_range) {
                            values(ib, r0_idx, itheta) = values(ib, r0_idx, theta0_idx);
                        }
                    });
        }
    }
};
