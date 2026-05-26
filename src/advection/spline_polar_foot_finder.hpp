// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "itimestepper.hpp"
#include "l_norm_tools.hpp"
#include "vector_index_tools.hpp"
#include "vector_mapper.hpp"

template <
        class GridR,
        class GridTheta,
        class X_adv,
        class Y_adv,
        class AdvDim1,
        class AdvDim2,
        class BSplinesR,
        class BSplinesTheta,
        class IdxRangeOperator,
        class SplineRThetaEvaluatorAdvection,
        class PseudoPhysicalToPhysicalMapping,
        class PseudoPhysicalToLogicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AdvecCoefField,
        class TimeStepper>
class ElementwiseSplinePolarFootFinder
{
    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using CoordRTheta = Coord<R, Theta>;

private:
    SplineRThetaEvaluatorAdvection m_evaluator_advection_field;

    PseudoPhysicalToPhysicalMapping m_pseudo_physical_to_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;

    TimeStepper m_time_stepper;

    AdvecCoefField m_advection_field_coefs;
    Coord<X_adv, Y_adv> m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;
    double m_dt;

public:
    ElementwiseSplinePolarFootFinder(
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            PseudoPhysicalToPhysicalMapping const& pseudo_physical_to_physical,
            PseudoPhysicalToLogicalMapping const& pseudo_physical_to_logical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            TimeStepper const& time_stepper,
            AdvecCoefField const& advection_field_coefs,
            Coord<X_adv, Y_adv> coord_centre,
            IdxRange<GridTheta> idx_range_theta,
            double dt)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_pseudo_physical_to_physical(pseudo_physical_to_physical)
        , m_pseudo_physical_to_logical(pseudo_physical_to_logical)
        , m_logical_to_pseudo_physical(logical_to_pseudo_physical)
        , m_advection_field_coefs(advection_field_coefs)
        , m_coord_centre(coord_centre)
        , m_idx_range_theta(idx_range_theta)
        , m_dt(dt)
    {
    }

    KOKKOS_FUNCTION CoordRTheta operator()(IdxOperator const idx) const
    {
        IdxBatch idx_batch(idx);
        IdxRTheta idx_rtheta(idx);
        // The function describing how the derivative of the evolve function is calculated.
        auto dy = [&](DVector<X_adv, Y_adv>& updated_advection_field, CoordRTheta const& foot) {
            DVector<AdvDim1, AdvDim2> updated_advection_field_adv_space;
            ddcHelper::get<AdvDim1>(updated_advection_field_adv_space)
                    = m_evaluator_advection_field(
                            foot,
                            get_const_field(
                                    ddcHelper::get<AdvDim1>(m_advection_field_coefs)[idx_batch]));
            ddcHelper::get<AdvDim2>(updated_advection_field_adv_space)
                    = m_evaluator_advection_field(
                            foot,
                            get_const_field(
                                    ddcHelper::get<AdvDim2>(m_advection_field_coefs)[idx_batch]));
            // Ensure coord is inside the domain as splines can't extrapolate
            // derivates (clamping)
            CoordRTheta advection_location_for_mapping(
                    Kokkos::min(ddc::select<R>(foot), ddc::discrete_space<BSplinesR>().rmax()),
                    ddc::select<Theta>(foot));
            updated_advection_field = to_vector_space<VectorIndexSet<X_adv, Y_adv>>(
                    m_pseudo_physical_to_physical,
                    advection_location_for_mapping,
                    updated_advection_field_adv_space);
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DVector<X_adv, Y_adv> const& advection_field,
                                   double dt) {
            Coord<X_adv, Y_adv> const coord_xy = m_logical_to_pseudo_physical(foot_rtheta);
            Coord<X_adv, Y_adv> const foot_xy = coord_xy - dt * advection_field;

            if (norm_inf(foot_xy - m_coord_centre) < 1e-15) {
                foot_rtheta = CoordRTheta(0, 0);
            } else {
                foot_rtheta = m_pseudo_physical_to_logical(foot_xy);
                ddc::select<Theta>(foot_rtheta) = ddcHelper::
                        restrict_to_idx_range(ddc::select<Theta>(foot_rtheta), m_idx_range_theta);
            }
        };

        CoordRTheta foot = ddc::coordinate(idx_rtheta);
        // Solve the characteristic equation
        m_time_stepper.update(foot, m_dt, dy, update_function);
        return foot;
    }
};

template <
        class GridR,
        class GridTheta,
        class X_adv,
        class Y_adv,
        class AdvDim1,
        class AdvDim2,
        class BSplinesR,
        class BSplinesTheta,
        class IdxRangeOperator,
        class SplineRThetaEvaluatorAdvection,
        class PseudoPhysicalToPhysicalMapping,
        class PseudoPhysicalToLogicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AdvecCoefFieldMem,
        class TimeStepper>
class ElementwiseSplinePolarFootFinderMem
{
    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using CoordRTheta = Coord<R, Theta>;
    using NonAlloc = ElementwiseSplinePolarFootFinder<
            GridR,
            GridTheta,
            X_adv,
            Y_adv,
            AdvDim1,
            AdvDim2,
            BSplinesR,
            BSplinesTheta,
            IdxRangeOperator,
            SplineRThetaEvaluatorAdvection,
            PseudoPhysicalToPhysicalMapping,
            PseudoPhysicalToLogicalMapping,
            LogicalToPseudoPhysicalMapping,
            typename AdvecCoefFieldMem::view_type,
            TimeStepper>;

private:
    SplineRThetaEvaluatorAdvection m_evaluator_advection_field;

    PseudoPhysicalToPhysicalMapping m_pseudo_physical_to_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;

    TimeStepper m_time_stepper;

    AdvecCoefFieldMem m_advection_field_coefs_alloc;
    Coord<X_adv, Y_adv> m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;

public:
    ElementwiseSplinePolarFootFinderMem(
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            PseudoPhysicalToPhysicalMapping const& pseudo_physical_to_physical,
            PseudoPhysicalToLogicalMapping const& pseudo_physical_to_logical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            TimeStepper const& time_stepper,
            AdvecCoefFieldMem&& advection_field_coefs,
            Coord<X_adv, Y_adv> coord_centre,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_pseudo_physical_to_physical(pseudo_physical_to_physical)
        , m_pseudo_physical_to_logical(pseudo_physical_to_logical)
        , m_logical_to_pseudo_physical(logical_to_pseudo_physical)
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_coord_centre(coord_centre)
        , m_idx_range_theta(idx_range_theta)
    {
    }

    NonAlloc operator()(double dt)
    {
        return NonAlloc(
                m_evaluator_advection_field,
                m_pseudo_physical_to_physical,
                m_pseudo_physical_to_logical,
                m_logical_to_pseudo_physical,
                m_time_stepper,
                get_const_field(m_advection_field_coefs_alloc),
                m_coord_centre,
                m_idx_range_theta,
                dt);
    }
};

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
        class IdxRangeBatched,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping,
        concepts::AnalyticalMapping LogicalToPseudoPhysicalMapping,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplinePolarFootFinder
{
    static_assert(is_timestepper_builder_v<TimeStepperBuilder>);
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

public:
    /// The continuous radial dimension.
    using GridR = typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1;
    /// The continuous poloidal dimension.
    using GridTheta = typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2;

    /// The type of the index range over which the operator works.
    using IdxRangeOperator = IdxRangeBatched;
    /// The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = typename SplineRThetaBuilderAdvection::memory_space;

    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    /// The continuous radial dimension.
    using VectorIndexSetAdvectionDims
            = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;

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

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
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

    using TimeStepper = typename TimeStepperBuilder::
            template time_stepper_t<CoordRTheta, DVector<X_adv, Y_adv>>;

    TimeStepperBuilder const& m_time_stepper_builder;

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

    using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
    using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;
    using ElementwiseOperator = ElementwiseSplinePolarFootFinderMem<
            GridR,
            GridTheta,
            X_adv,
            Y_adv,
            AdvDim1,
            AdvDim2,
            BSplinesR,
            BSplinesTheta,
            IdxRangeOperator,
            SplineRThetaEvaluatorAdvection,
            PseudoPhysicalToPhysicalMapping,
            PseudoPhysicalToLogicalMapping,
            LogicalToPseudoPhysicalMapping,
            DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>,
            TimeStepper>;

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
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical_mapping,
            SplineRThetaBuilderAdvection const& builder_advection_field,
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            double epsilon = 1e-12)
        : m_time_stepper_builder(time_stepper_builder)
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
    ElementwiseOperator operator()(
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field)
    {
        static_assert(ddc::type_seq_size_v<VectorIndexSetAdvectionDims> == 2);

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

        ElementwiseOperator elementwise(
                m_evaluator_advection_field,
                m_pseudo_physical_to_physical,
                m_pseudo_physical_to_logical,
                m_logical_to_pseudo_physical,
                time_stepper,
                std::move(advection_field_coefs),
                m_logical_to_pseudo_physical(CoordRTheta(0, 0)),
                idx_range_theta);
        return elementwise;
    }
    void operator()(
            CFieldFeet feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt)
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
                m_pseudo_physical_to_physical,
                m_pseudo_physical_to_logical,
                m_logical_to_pseudo_physical,
                time_stepper,
                std::move(advection_field_coefs),
                m_logical_to_pseudo_physical(CoordRTheta(0, 0)),
                idx_range_theta);

        auto elementwise = elementwise_mem(dt);

        // Compute the characteristic feet at t^n:
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) { feet(idx) = elementwise(idx); });

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
