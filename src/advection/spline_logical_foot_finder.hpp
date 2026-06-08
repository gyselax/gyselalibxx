// SPDX-License-Identifier: MIT
#pragma once
#include <source_location>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "itimestepper.hpp"
#include "l_norm_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

/**
 * @brief A device-callable functor that finds the foot of the characteristic at a single
 *        grid point in polar coordinates using spline interpolation.
 *
 * This class holds non-owning views of all the objects needed to solve the
 * characteristic equation at a single index, and so it can be copied to and
 * called from the device. It is the innermost computation unit of
 * @ref SplineLogicalFootFinder.
 *
 * Unlike @ref ElementwiseSplinePolarFootFinder, the characteristic equation is integrated
 * directly in the logical @f$(r, \theta)@f$ coordinate system. No coordinate transformation
 * to a pseudo-Cartesian domain is performed. This means the advection field must be provided
 * in the contravariant logical basis @f$(A^r, A^\theta)@f$.
 *
 * If the characteristic foot has @f$ r < 0 @f$ after a time step, it is reflected
 * through the O-point: @f$ r \to -r @f$, @f$ \theta \to \theta + \pi @f$.
 *
 * @tparam GridR
 *      The discrete radial dimension.
 * @tparam GridTheta
 *      The discrete poloidal dimension.
 * @tparam BSplinesR
 *      The B-spline basis in the radial direction.
 * @tparam BSplinesTheta
 *      The B-spline basis in the poloidal direction.
 * @tparam IdxRangeOperator
 *      The full index range over which the operator acts (may include batch dimensions).
 * @tparam SplineRThetaEvaluatorAdvection
 *      The evaluator used to evaluate the spline representation of the advection field.
 * @tparam AdvecCoefField
 *      A non-owning field (view) type holding the spline coefficients of the advection field.
 * @tparam TimeStepper
 *      The time integration method used to solve the characteristic equation.
 *
 * @see ElementwiseSplineLogicalFootFinderMem
 * @see SplineLogicalFootFinder
 */
template <
        class GridR,
        class GridTheta,
        class BSplinesR,
        class BSplinesTheta,
        class IdxRangeOperator,
        class SplineRThetaEvaluatorAdvection,
        class AdvecCoefField,
        class TimeStepper>
class ElementwiseSplineLogicalFootFinder
{
    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    /// First advection dimension (radial, in logical basis).
    using AdvDim1 = R;
    /// Second advection dimension (poloidal, in logical basis).
    using AdvDim2 = Theta;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using CoordRTheta = Coord<R, Theta>;

private:
    SplineRThetaEvaluatorAdvection m_evaluator_advection_field;
    TimeStepper m_time_stepper;
    AdvecCoefField m_advection_field_coefs;
    IdxRange<GridTheta> m_idx_range_theta;
    double m_dt;

public:
    /**
     * @brief Construct an ElementwiseSplineLogicalFootFinder.
     *
     * @param[in] evaluator_advection_field
     *      The evaluator for the spline representation of the advection field.
     * @param[in] time_stepper
     *      The time integration method.
     * @param[in] advection_field_coefs
     *      A non-owning view of the pre-built spline coefficients of the advection field.
     * @param[in] idx_range_theta
     *      The poloidal index range, used to wrap the angular coordinate into
     *      the periodic domain after each time step.
     * @param[in] dt
     *      The time step for the characteristic equation.
     */
    ElementwiseSplineLogicalFootFinder(
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            TimeStepper const& time_stepper,
            AdvecCoefField const& advection_field_coefs,
            IdxRange<GridTheta> idx_range_theta,
            double dt)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_time_stepper(time_stepper)
        , m_advection_field_coefs(advection_field_coefs)
        , m_idx_range_theta(idx_range_theta)
        , m_dt(dt)
    {
    }

    /**
     * @brief Find the foot of the characteristic at a single grid index.
     *
     * Solves the characteristic equation over @f$ dt @f$ using the stored
     * time stepper and returns the resulting @f$ (r, \theta) @f$ coordinate.
     * The integration is performed directly in the logical coordinate system.
     *
     * @param[in] idx
     *      The operator index, encoding both batch and @f$ (r, \theta) @f$ indices.
     *
     * @return The @f$ (r, \theta) @f$ coordinate of the characteristic foot.
     */
    KOKKOS_FUNCTION CoordRTheta operator()(IdxOperator const idx) const
    {
        IdxBatch idx_batch(idx);
        IdxRTheta idx_rtheta(idx);

        // Evaluate the advection field in logical coordinates at the current foot.
        auto dy = [&](DVector<R, Theta>& adv, CoordRTheta const& foot) {
            ddcHelper::get<AdvDim1>(adv) = m_evaluator_advection_field(
                    foot,
                    ddcHelper::get<AdvDim1>(m_advection_field_coefs)[idx_batch]);
            ddcHelper::get<AdvDim2>(adv) = m_evaluator_advection_field(
                    foot,
                    ddcHelper::get<AdvDim2>(m_advection_field_coefs)[idx_batch]);
        };

        // Update the foot directly in logical coordinates.
        auto update_function = [&](CoordRTheta& foot, DVector<R, Theta> const& adv, double dt) {
            foot -= dt * adv;

            // O-point reflection: if r goes negative, reflect through the origin.
            double foot_r = ddc::get<R>(foot);
            double foot_theta = ddc::get<Theta>(foot);
            int negative_reflexion = static_cast<int>(foot_r < 0);
            foot = CoordRTheta(
                    (1 - 2 * negative_reflexion) * foot_r,
                    foot_theta + M_PI * negative_reflexion);
            // Wrap theta into the periodic domain.
            ddc::select<Theta>(foot)
                    = ddcHelper::restrict_to_idx_range(ddc::select<Theta>(foot), m_idx_range_theta);
        };

        CoordRTheta foot = ddc::coordinate(idx_rtheta);
        m_time_stepper.update(foot, m_dt, dy, update_function);
        return foot;
    }
};


/**
 * @brief The owning counterpart to ElementwiseSplineLogicalFootFinder.
 *
 * Allocates and stores the spline coefficients of the advection field on the
 * appropriate memory space. Calling @c operator()(dt) returns a non-owning
 * @ref ElementwiseSplineLogicalFootFinder configured for the given time step, which
 * can then be copied to and called from the device.
 *
 * @tparam GridR
 *      The discrete radial dimension.
 * @tparam GridTheta
 *      The discrete poloidal dimension.
 * @tparam BSplinesR
 *      The B-spline basis in the radial direction.
 * @tparam BSplinesTheta
 *      The B-spline basis in the poloidal direction.
 * @tparam IdxRangeOperator
 *      The full index range over which the operator acts (may include batch dimensions).
 * @tparam SplineRThetaEvaluatorAdvection
 *      The evaluator used to evaluate the spline representation of the advection field.
 * @tparam AdvecCoefFieldMem
 *      An owning field type holding the spline coefficients of the advection field.
 * @tparam TimeStepper
 *      The time integration method used to solve the characteristic equation.
 *
 * @see ElementwiseSplineLogicalFootFinder
 * @see SplineLogicalFootFinder
 */
template <
        class GridR,
        class GridTheta,
        class BSplinesR,
        class BSplinesTheta,
        class IdxRangeOperator,
        class SplineRThetaEvaluatorAdvection,
        class AdvecCoefFieldMem,
        class TimeStepper>
class ElementwiseSplineLogicalFootFinderMem
{
    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;

public:
    /// The non-owning operator that can be used on GPU.
    using GPUCompat = ElementwiseSplineLogicalFootFinder<
            GridR,
            GridTheta,
            BSplinesR,
            BSplinesTheta,
            IdxRangeOperator,
            SplineRThetaEvaluatorAdvection,
            typename AdvecCoefFieldMem::view_type,
            TimeStepper>;

private:
    SplineRThetaEvaluatorAdvection m_evaluator_advection_field;
    TimeStepper m_time_stepper;
    AdvecCoefFieldMem m_advection_field_coefs_alloc;
    IdxRange<GridTheta> m_idx_range_theta;

public:
    /**
     * @brief Construct an ElementwiseSplineLogicalFootFinderMem.
     *
     * @param[in] evaluator_advection_field
     *      The evaluator for the spline representation of the advection field.
     * @param[in] time_stepper
     *      The time integration method.
     * @param[in] advection_field_coefs
     *      The spline coefficients of the advection field. Ownership is transferred in.
     * @param[in] idx_range_theta
     *      The poloidal index range, used to wrap the angular coordinate into
     *      the periodic domain after each time step.
     */
    ElementwiseSplineLogicalFootFinderMem(
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            TimeStepper const& time_stepper,
            AdvecCoefFieldMem&& advection_field_coefs,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_time_stepper(time_stepper)
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_idx_range_theta(idx_range_theta)
    {
    }

    /**
     * @brief Create an ElementwiseSplineLogicalFootFinder for the given time step.
     *
     * Returns a non-owning @ref ElementwiseSplineLogicalFootFinder that holds views of
     * the stored spline coefficients and is configured for time step @f$ dt @f$.
     * The returned object can be copied to and called from the device.
     *
     * @param[in] dt
     *      The time step for the characteristic equation.
     *
     * @return A view-based ElementwiseSplineLogicalFootFinder for the given time step.
     */
    GPUCompat operator()(double dt)
    {
        return GPUCompat(
                m_evaluator_advection_field,
                m_time_stepper,
                get_const_field(m_advection_field_coefs_alloc),
                m_idx_range_theta,
                dt);
    }
};


/**
 * @brief A class to find the foot of the characteristics on the @f$ (r,\theta) @f$
 * plane, integrating directly in the logical coordinate system.
 *
 * Unlike @ref SplinePolarFootFinder, this class does not convert to a pseudo-Cartesian
 * domain. The characteristic equation is solved directly in the logical @f$ (r,\theta) @f$
 * coordinate system. As a consequence, the advection field must be provided in the
 * contravariant logical basis @f$ (A^r, A^\theta) @f$.
 *
 * If the characteristic foot has @f$ r < 0 @f$ after a time step, it is reflected
 * through the O-point: @f$ r \to -r @f$, @f$ \theta \to \theta + \pi @f$.
 *
 * @tparam IdxRangeBatched
 *      The full index range over which the operator acts (may include batch dimensions).
 * @tparam TimeStepperBuilder
 *      A time stepper builder indicating which time integration method should be
 *      applied to solve the characteristic equation.
 * @tparam SplineRThetaBuilderAdvection
 *      A 2D SplineBuilder to construct a spline on a polar domain.
 * @tparam SplineRThetaEvaluatorAdvection
 *      A 2D SplineEvaluator to evaluate a spline on a polar domain.
 *      A boundary condition must be provided in case the foot of the characteristic
 *      is found outside the domain.
 *
 * @see BslAdvectionPolar
 * @see SplinePolarFootFinder
 */
template <
        class IdxRangeBatched,
        class TimeStepperBuilder,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplineLogicalFootFinder
{
    static_assert(is_timestepper_builder_v<TimeStepperBuilder>);
    static_assert(std::is_same_v<
                  typename SplineRThetaBuilderAdvection::memory_space,
                  typename SplineRThetaEvaluatorAdvection::memory_space>);
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
    /// The discrete radial dimension.
    using GridR = typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1;
    /// The discrete poloidal dimension.
    using GridTheta = typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2;

    /// The type of the index range over which the operator works.
    using IdxRangeOperator = IdxRangeBatched;
    /// The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = typename SplineRThetaBuilderAdvection::memory_space;

    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    /// The logical coordinate type.
    using CoordRTheta = Coord<R, Theta>;

    /**
     * @brief The dimensions of the advection field: the logical contravariant basis
     * @f$ (R, \Theta) @f$.
     */
    using VectorIndexSetAdvectionDims = ddc::to_type_seq_t<CoordRTheta>;

    /// @brief Execution space.
    using ExecSpace = typename SplineRThetaBuilderAdvection::exec_space;

private:
    using BSplinesR = typename SplineRThetaBuilderAdvection::bsplines_type1;
    using BSplinesTheta = typename SplineRThetaBuilderAdvection::bsplines_type2;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;

    using IdxRangeSplineBatched
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<IdxRangeOperator>,
                    ddc::detail::TypeSeq<GridR, GridTheta>,
                    ddc::detail::TypeSeq<BSplinesR, BSplinesTheta>>>;

    using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
    using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

    using TimeStepper =
            typename TimeStepperBuilder::template time_stepper_t<CoordRTheta, DVector<R, Theta>>;

    TimeStepperBuilder const& m_time_stepper_builder;
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
     * @brief The type of a vector field defined on the logical basis at every
     * grid point, saved on a compatible memory space.
     */
    using DVectorFieldAdvection
            = DVectorField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>;

    /**
     * @brief The type of a constant vector field defined on the logical basis
     * at every grid point, saved on a compatible memory space.
     */
    using DVectorConstFieldAdvection
            = DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>;

    /// The operator returned by operator() which calculates the feet elementwise.
    using ElementwiseOperator = ElementwiseSplineLogicalFootFinderMem<
            GridR,
            GridTheta,
            BSplinesR,
            BSplinesTheta,
            IdxRangeOperator,
            SplineRThetaEvaluatorAdvection,
            DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>,
            TimeStepper>;

    /**
     * @brief Instantiate a foot finder that integrates directly in the logical domain.
     *
     * @param[in] idx_range_operator
     *      The index range on which the operator should act. Used only for
     *      template argument deduction of @p IdxRangeBatched.
     * @param[in] time_stepper_builder
     *      A builder for the time integration method used for the
     *      characteristic equation.
     * @param[in] builder_advection_field
     *      The spline builder which computes the spline representation
     *      of the advection field.
     * @param[in] evaluator_advection_field
     *      The B-splines evaluator to evaluate the advection field.
     */
    SplineLogicalFootFinder(
            IdxRangeBatched const& idx_range_operator,
            TimeStepperBuilder const& time_stepper_builder,
            SplineRThetaBuilderAdvection const& builder_advection_field,
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field)
        : m_time_stepper_builder(time_stepper_builder)
        , m_builder_advection_field(builder_advection_field)
        , m_evaluator_advection_field(evaluator_advection_field)
    {
    }

    /**
     * @brief Get an elementwise operator providing a GPU copyable functor capable of
     * calculating the feet of the characteristics.
     *
     * Computes the B-spline coefficients of the advection field and returns a
     * device-callable @ref ElementwiseSplineLogicalFootFinderMem.
     *
     * @param[in] advection_field
     *      The advection field in the logical contravariant basis.
     *
     * @returns An elementwise operator providing a GPU copyable functor capable of
     *      calculating the feet of the characteristics.
     */
    ElementwiseOperator operator()(
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field) const
    {
        DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>
                advection_field_coefs(m_builder_advection_field.batched_spline_domain(
                        get_idx_range(advection_field)));

        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));

        IdxRangeTheta idx_range_theta(get_idx_range(advection_field));
        TimeStepper time_stepper = m_time_stepper_builder.template preallocate<TimeStepper>();

        return ElementwiseOperator(
                m_evaluator_advection_field,
                time_stepper,
                std::move(advection_field_coefs),
                idx_range_theta);
    }

    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * Computes the B-spline coefficients of the advection field, then integrates
     * the characteristic equation for each grid point.
     *
     * @param[in, out] feet
     *      On input: the mesh points.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the logical contravariant basis.
     * @param[in] dt
     *      The time step.
     */
    void operator()(
            CFieldFeet feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt) const
    {
        DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>
                advection_field_coefs(m_builder_advection_field.batched_spline_domain(
                        get_idx_range(advection_field)));

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
                time_stepper,
                std::move(advection_field_coefs),
                idx_range_theta);

        typename ElementwiseOperator::GPUCompat elementwise = elementwise_mem(dt);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) { feet(idx) = elementwise(idx); });

        unify_value_at_centre_pt(feet);
        is_unified(feet);
    }

    /**
     * @brief Check if the values at the centre point are the same.
     *
     *  For polar geometry, to ensure continuity at the centre point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  This function checks if for @f$ r= 0 @f$, the values @f$ \forall \theta @f$ are the same.
     *
     *  @param[in] values
     *      A table of values we want to check if the centre point has
     *      a unique value.
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
     *  As the computation of the values of a table can introduce machine errors,
     *  this function is useful to reset the values at the central point to
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
