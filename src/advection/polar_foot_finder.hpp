// SPDX-License-Identifier: MIT
#pragma once
#include <source_location>

#include <ddc/ddc.hpp>

#include "polar_foot_finders/elementwise_choice.hpp"
#include "polar_foot_finders/logical_advection_logical_foot_finder.hpp"
#include "polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp"
#include "polar_foot_finders/physical_advection_physical_foot_finder.hpp"
#include "polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "i_interpolation.hpp"
#include "l_norm_tools.hpp"
#include "vector_index_tools.hpp"

/**
 * @brief Operator for finding the feet of the characteristics on a polar slice.
 *
 * Calculates the interpolation representation of the advection field and uses it together
 * with a time-stepping method to solve the characteristic equation. The space in
 * which the advection field is expressed and the space in which foot-finding is
 * performed are selected at compile time via @p FFSpace and @p AFSpace.
 *
 * @tparam FFSpace     The space in which the foot of the characteristic is found.
 * @tparam AFSpace     The space in which the advection field is expressed.
 * @tparam LogicalToPhysicalMapping
 *      The mapping from the logical @f$ (r, \theta) @f$ domain to the physical
 *      @f$ (x, y) @f$ domain.
 * @tparam IdxRangeBatched
 *      The batched index range over which the operator works.
 * @tparam TimeStepperBuilder
 *      The factory type for the time integration method.
 * @tparam RThetaAdvectionInterpolator
 *      The interpolator for the advection field.
 */
template <
        FootFindingSpace FFSpace,
        AdvectionFieldSpace AFSpace,
        concepts::Mapping LogicalToPhysicalMapping,
        class IdxRangeBatched,
        class TimeStepperBuilder,
        concepts::Interpolation RThetaAdvectionInterpolator>
class PolarFootFinder
{
    static_assert(
            FFSpace != FootFindingSpace::PHYSICAL
                    || is_analytical_mapping_v<LogicalToPhysicalMapping>,
            "It is not possible to find the foot of the characteristic in Physical space as there "
            "is no way to return to Logical space once the foot is calculated.");

private:
    using RThetaAdvectionBuilder = typename RThetaAdvectionInterpolator::BuilderType;
    using RThetaAdvectionEvaluator = typename RThetaAdvectionInterpolator::EvaluatorType;

public:
    /// The continuous radial dimension.
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    /// The continuous poloidal dimension.
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    /// The memory space where fields are stored (e.g. @c Kokkos::HostSpace or a GPU space).
    using memory_space = typename RThetaAdvectionBuilder::memory_space;
    /// The execution space where kernels are launched.
    using ExecSpace = typename RThetaAdvectionBuilder::exec_space;

private:
    using LogicalSpace = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordArg>;
    using PhysicalSpace = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;

    using X = ddc::type_seq_element_t<0, PhysicalSpace>;
    using Y = ddc::type_seq_element_t<1, PhysicalSpace>;

public:
    /// The first dimension on which the advection field is defined.
    using AdvDim1 = std::conditional_t<AFSpace == AdvectionFieldSpace::PHYSICAL, X, R>;
    /// The second dimension on which the advection field is defined.
    using AdvDim2 = std::conditional_t<AFSpace == AdvectionFieldSpace::PHYSICAL, Y, Theta>;

private:
    using PolarGrid = ddc::to_type_seq_t<typename InterpolationBuilderTraits<
            RThetaAdvectionBuilder>::interpolation_idx_range_type>;
    using GridR = find_grid_t<R, PolarGrid>;
    using GridTheta = find_grid_t<Theta, PolarGrid>;

    using IdxRangeCoeffBatched = typename InterpolationBuilderTraits<
            RThetaAdvectionBuilder>::batched_basis_idx_range_type<IdxRangeBatched>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatched, GridR, GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxOperator = typename IdxRangeBatched::discrete_element_type;

    using CoordRTheta = Coord<R, Theta>;

    using AdvecCoefField
            = DVectorFieldMem<IdxRangeCoeffBatched, VectorIndexSet<AdvDim1, AdvDim2>, memory_space>;

public:
    /// The operator returned by operator() which calculates the feet elementwise.
    using ElementwiseOperator = typename polar_foot_finder_details::ElementwiseChoice<
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
    /// The type of the index range over which the operator works.
    using IdxRangeOperator = IdxRangeBatched;
    /// The type of the memory space where the field is saved (CPU vs GPU).
    /// A mutable field of @f$ (r, \theta) @f$ coordinates over the batched operator domain.
    using CFieldFeet = Field<CoordRTheta, IdxRangeBatched, memory_space>;
    /// A read-only field of @f$ (r, \theta) @f$ coordinates over the batched operator domain.
    using CConstFieldFeet = ConstField<CoordRTheta, IdxRangeBatched, memory_space>;

private:
    TimeStepperBuilder const& m_time_stepper_builder;
    LogicalToPhysicalMapping m_logical_to_physical;
    RThetaAdvectionBuilder const& m_builder_advection_field;
    RThetaAdvectionEvaluator const& m_evaluator_advection_field;
    Coord<X_pC, Y_pC> m_coord_centre_pc;
    double m_epsilon;

public:
    /**
     * @brief Construct a PolarFootFinder.
     *
     * @param[in] time_stepper_builder
     *      A builder for the time integration method.
     * @param[in] logical_to_physical
     *      The mapping from the logical domain to the physical domain.
     * @param[in] interpolator_advection_field
     *      An interpolator to build and evaluate an approximation of the
     *      advection field.
     * @param[in] coord_centre_pc
     *      The coordinate of the polar centre in the pseudo-Cartesian domain
     *      @f$ (X_{pC}, Y_{pC}) @f$. Ignored for @c LOGICAL foot finding.
     * @param[in] epsilon
     *      @f$ \varepsilon @f$ parameter for the CombinedMapping linearisation
     *      near the O-point. Only used for @c PHYSICAL advection field.
     */
    PolarFootFinder(
            TimeStepperBuilder const& time_stepper_builder,
            LogicalToPhysicalMapping const& logical_to_physical,
            RThetaAdvectionInterpolator const& interpolator_advection_field,
            Coord<X_pC, Y_pC> coord_centre_pc = Coord<X_pC, Y_pC>(0, 0),
            double epsilon = 1e-12)
        : m_time_stepper_builder(time_stepper_builder)
        , m_logical_to_physical(logical_to_physical)
        , m_builder_advection_field(interpolator_advection_field.get_builder())
        , m_evaluator_advection_field(interpolator_advection_field.get_evaluator())
        , m_coord_centre_pc(coord_centre_pc)
        , m_epsilon(epsilon)
    {
    }

    /**
     * @brief Get an elementwise operator capable of calculating the feet of the characteristics.
     *
     * Computes the interpolation coefficients of the advection field, then packages them together
     * with the mappings and time stepper into an @ref ElementwiseOperator. Calling
     * @c operator()(dt) on the returned object yields a GPU-copyable functor.
     *
     * @param[in] advection_field
     *      The advection field in the domain selected by @c AFSpace.
     * @returns An @ref ElementwiseOperator for the given advection field.
     */
    ElementwiseOperator operator()(
            DVectorConstField<IdxRangeBatched, VectorIndexSet<AdvDim1, AdvDim2>, memory_space>
                    advection_field) const
    {
        AdvecCoefField advection_field_coefs(
                batched_basis_idx_range(m_builder_advection_field, get_idx_range(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));
        IdxRangeTheta idx_range_theta(get_idx_range(advection_field));
        return build_elementwise(std::move(advection_field_coefs), idx_range_theta);
    }

    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * Computes the coefficients of the advection field, solves the characteristic
     * equation over @f$ dt @f$ at every grid point, and writes the resulting feet in-place.
     *
     * @param[in, out] feet
     *      On input: the mesh points. On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the domain selected by @c AFSpace.
     * @param[in] dt
     *      The time step.
     */
    void operator()(
            CFieldFeet feet,
            DVectorConstField<IdxRangeBatched, VectorIndexSet<AdvDim1, AdvDim2>, memory_space>
                    advection_field,
            double dt) const
    {
        AdvecCoefField advection_field_coefs(
                batched_basis_idx_range(m_builder_advection_field, get_idx_range(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));
        IdxRangeTheta idx_range_theta(get_idx_range(advection_field));

        ElementwiseOperator elementwise_mem
                = build_elementwise(std::move(advection_field_coefs), idx_range_theta);
        typename ElementwiseOperator::GPUCompat elementwise = elementwise_mem(dt);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) { feet(idx) = elementwise(idx); });
    }

private:
    ElementwiseOperator build_elementwise(
            AdvecCoefField&& advection_field_coefs,
            IdxRangeTheta idx_range_theta) const
    {
        if constexpr (
                FFSpace == FootFindingSpace::PSEUDO_PHYSICAL
                && AFSpace == AdvectionFieldSpace::PHYSICAL) {
            // PSEUDO_PHYSICAL/PHYSICAL: uses global pseudo-Cartesian (X_pC, Y_pC) space.
            // The CombinedMapping inside the Mem class requires epsilon.
            return ElementwiseOperator(
                    m_evaluator_advection_field,
                    m_logical_to_physical,
                    m_time_stepper_builder,
                    std::move(advection_field_coefs),
                    m_coord_centre_pc,
                    idx_range_theta,
                    m_epsilon);
        } else if constexpr (FFSpace == FootFindingSpace::PHYSICAL) {
            // PHYSICAL/LOGICAL: foot finding in physical (X, Y) space.
            // coord_centre type is Coord<X, Y>, computed from the physical mapping.
            return ElementwiseOperator(
                    m_evaluator_advection_field,
                    m_logical_to_physical,
                    m_time_stepper_builder,
                    std::move(advection_field_coefs),
                    m_logical_to_physical(CoordRTheta(0, 0)),
                    idx_range_theta);
        } else {
            // LOGICAL/LOGICAL or PSEUDO_PHYSICAL/LOGICAL:
            // coord_centre is Coord<X_pC, Y_pC> (unused for LOGICAL/LOGICAL).
            return ElementwiseOperator(
                    m_evaluator_advection_field,
                    m_logical_to_physical,
                    m_time_stepper_builder,
                    std::move(advection_field_coefs),
                    m_coord_centre_pc,
                    idx_range_theta);
        }
    }
};

/**
 * @brief Construct a PolarFootFinder, deducing all type template parameters from the arguments.
 *
 * @c FFSpace and @c AFSpace must still be provided explicitly (they are enum values, not types).
 * All remaining template parameters are deduced from the function arguments, eliminating the
 * need for @c decltype at the call site.
 *
 * @param[in] time_stepper  A builder for the time integration method.
 * @param[in] mapping       The mapping from the logical domain to the physical domain.
 * @param[in] idx_range     The batched index range over which the operator works (used only for
 *                          type deduction; its value is not forwarded to the constructor).
 * @param[in] interpolator_advection_field
 *                          An interpolator to build and evaluate an approximation of the
 *                          advection field.
 * @param[in] coord_centre  The polar-centre coordinate in pseudo-Cartesian space.
 * @param[in] epsilon       Linearisation parameter near the O-point.
 */
template <
        FootFindingSpace FFSpace,
        AdvectionFieldSpace AFSpace,
        concepts::Mapping LogicalToPhysicalMapping,
        class IdxRangeBatched,
        class TimeStepperBuilder,
        concepts::Interpolation RThetaAdvectionInterpolator>
auto make_polar_foot_finder(
        TimeStepperBuilder const& time_stepper,
        LogicalToPhysicalMapping const& mapping,
        [[maybe_unused]] IdxRangeBatched const& idx_range,
        RThetaAdvectionInterpolator const& interpolator_advection_field,
        Coord<X_pC, Y_pC> coord_centre = Coord<X_pC, Y_pC>(0, 0),
        double epsilon = 1e-12)
{
    return PolarFootFinder<
            FFSpace,
            AFSpace,
            LogicalToPhysicalMapping,
            IdxRangeBatched,
            TimeStepperBuilder,
            RThetaAdvectionInterpolator>(
            time_stepper,
            mapping,
            interpolator_advection_field,
            coord_centre,
            epsilon);
}
