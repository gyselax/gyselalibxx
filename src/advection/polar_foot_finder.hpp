#pragma once
#include "polar_foot_finders/elementwise_choice.hpp"

template <
        FootFindingSpace FFSpace,
        AdvectionFieldSpace AFSpace,
        concepts::Mapping LogicalToPhysicalMapping,
        class IdxRangeBatched,
        class TimeStepperBuilder,
        class RThetaAdvectionBuilder,
        class RThetaAdvectionEvaluator>
class PolarFootFinder
{
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    using LogicalSpace = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordArg>;
    using PhysicalSpace = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;

    using X = ddc::type_seq_element_t<0, PhysicalSpace>;
    using Y = ddc::type_seq_element_t<1, PhysicalSpace>;

    using PolarGrid
            = ddc::to_type_seq_t<typename RThetaAdvectionBuilder::interpolation_domain_type>;
    using GridR = find_grid_t<R, PolarGrid>;
    using GridTheta = find_grid_t<Theta, PolarGrid>;

    using AdvDim1 = std::conditional_t<AFSpace == PHYSICAL, X, R>;
    using AdvDim2 = std::conditional_t<AFSpace == PHYSICAL, Y, Theta>;

    //using FFDim1 = std::conditional_t<FFSpace == PHYSICAL, X, std::conditional_t<FFSpace == PSEUDO_PHYSICAL, X_pC, R>>;
    //using FFDim2 = std::conditional_t<FFSpace == PHYSICAL, Y, std::conditional_t<FFSpace == PSEUDO_PHYSICAL, Y_pC, Theta>>;

    static_assert(FFSpace != PHYSICAL || is_analytical_mapping_v<LogicalToPhysicalMapping>);

    /// The operator returned by operator() which calculates the feet elementwise.
    using ElementwiseOperator = ElementwiseChoice<
            FFSpace,
            AFSpace,
            GridR,
            GridTheta,
            IdxRangeBatched,
            RThetaAdvectionEvaluator,
            AdvecCoefField,
            TimeStepperBuilder,
            LogicalToPhysicalMapping>::type;

public:
    /**
     * @brief Get an elementwise operator providing a GPU copyable functor capable of
     * calculating the feet of the characteristics.
     *
     * From the advection field in the physical domain, compute the advection field
     * in the right domain an compute its B-splines coefficients.
     * Then, use the given time integration method (time_stepper) to solve the
     * characteristic equation over @f$ dt @f$.
     *
     * @param[in] advection_field
     *      The advection field in the chosen domain.
     *
     * @returns An elementwise operator providing a GPU copyable functor capable of
     *      calculating the feet of the characteristics.
     */
    ElementwiseOperator operator()(
            DVectorConstField<IdxRangeOperator, VectorIndexSet<AdvDim1, AdvDim2>, memory_space>
                    advection_field) const
    {
        static_assert(ddc::type_seq_size_v<VectorIndexSetAdvectionDims> == 2);

        DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSet<AdvDim1, AdvDim2>, memory_space>
                advection_field_coefs(m_builder_advection_field.batched_spline_domain(
                        get_idx_range(advection_field)));

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));

        IdxRangeTheta idx_range_theta(get_idx_range(advection_field));

        ElementwiseOperator elementwise(
                m_evaluator_advection_field,
                m_pseudo_physical_to_advection,
                m_pseudo_physical_to_logical,
                m_logical_to_pseudo_physical,
                m_time_stepper_builder,
                std::move(advection_field_coefs),
                m_logical_to_pseudo_physical(CoordRTheta(0, 0)),
                idx_range_theta);
        return elementwise;
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
            double dt) const
    {
        static_assert(ddc::type_seq_size_v<VectorIndexSetAdvectionDims> == 2);
        using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
        using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

        DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>
                advection_field_coefs(m_builder_advection_field.batched_spline_domain(
                        get_idx_range(advection_field)));

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));

        IdxRangeTheta idx_range_theta(get_idx_range(advection_field));

        TimeStepper time_stepper = m_time_stepper_builder.template preallocate<TimeStepper>();

        ElementwiseOperator elementwise_mem(
                m_evaluator_advection_field,
                m_pseudo_physical_to_advection,
                m_pseudo_physical_to_logical,
                m_logical_to_pseudo_physical,
                time_stepper,
                std::move(advection_field_coefs),
                m_logical_to_pseudo_physical(CoordRTheta(0, 0)),
                idx_range_theta);

        typename ElementwiseOperator::GPUCompat elementwise = elementwise_mem(dt);

        // Compute the characteristic feet at t^n:
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) { feet(idx) = elementwise(idx); });

        //// Treatment to conserve the C0 property of the advected function:
        //unify_value_at_centre_pt(feet);
        //// Test if the values are the same at the centre point
        //is_unified(feet);
    }
};
