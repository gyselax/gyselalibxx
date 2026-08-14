// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "l_norm_tools.hpp"
#include "tensor.hpp"
#include "type_seq_tools.hpp"
#include "vector_field.hpp"

namespace polar_foot_finder_details {

/**
 * @brief GPU-callable functor that finds characteristic feet for logical-space
 *        advection with foot-finding in pseudo-physical space.
 *
 * The advection field is expressed in logical @f$ (r, \theta) @f$ coordinates and
 * is converted to pseudo-physical @f$ (X_{pC}, Y_{pC}) @f$ space before the
 * characteristic equation is integrated. Constructed by
 * @ref ElementwiseLogicalAdvPseudoPhysFootFinderMem.
 */
template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepper,
        class LogicalToPseudoPhysicalMapping>
class ElementwiseLogicalAdvPseudoPhysFootFinder
{
    using PseudoPhysicalToLogicalMapping = inverse_mapping_t<LogicalToPseudoPhysicalMapping>;
    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;
    using BSplinesR = find_grid_t<
            R,
            ddc::to_type_seq_t<typename RThetaAdvectionEvaluator::spline_domain_type>>;
    using IdxRTheta = Idx<GridR, GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeOperator, GridR, GridTheta>;
    using IdxOperator = typename IdxRangeOperator::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using CoordRTheta = Coord<R, Theta>;
    using CoordX_pcY_pc = typename PseudoPhysicalToLogicalMapping::CoordArg;
    using PseudoCartBasis = ddc::to_type_seq_t<CoordX_pcY_pc>;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;

    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;

    TimeStepper m_time_stepper;

    AdvecCoefField m_advection_field_coefs;
    CoordX_pcY_pc m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;
    double m_dt;

public:
    /**
     * @brief Construct an ElementwiseLogicalAdvPseudoPhysFootFinder.
     *
     * @param[in] evaluator_advection_field
     *      The evaluator for the spline representation of the advection field.
     * @param[in] pseudo_physical_to_logical
     *      The mapping from the pseudo-physical domain to the logical domain.
     * @param[in] logical_to_pseudo_physical
     *      The mapping from the logical domain to the pseudo-physical domain.
     * @param[in] time_stepper
     *      The time integration method.
     * @param[in] advection_field_coefs
     *      A view of the spline coefficients of the advection field.
     * @param[in] coord_centre
     *      The coordinate of the polar centre in the pseudo-physical domain,
     *      used to handle the degenerate point at @f$ r = 0 @f$.
     * @param[in] idx_range_theta
     *      The poloidal index range, used to wrap the angular coordinate into
     *      the periodic domain after each time step.
     * @param[in] dt
     *      The time step for the characteristic equation.
     */
    ElementwiseLogicalAdvPseudoPhysFootFinder(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            PseudoPhysicalToLogicalMapping const& pseudo_physical_to_logical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            TimeStepper const& time_stepper,
            AdvecCoefField const& advection_field_coefs,
            CoordX_pcY_pc coord_centre,
            IdxRange<GridTheta> idx_range_theta,
            double dt)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_pseudo_physical_to_logical(pseudo_physical_to_logical)
        , m_logical_to_pseudo_physical(logical_to_pseudo_physical)
        , m_time_stepper(time_stepper)
        , m_advection_field_coefs(advection_field_coefs)
        , m_coord_centre(coord_centre)
        , m_idx_range_theta(idx_range_theta)
        , m_dt(dt)
    {
    }

    /**
     * @brief Find the foot of the characteristic at a single grid index.
     *
     * Solves the characteristic equation over @f$ dt @f$ using the stored
     * time stepper and returns the resulting @f$ (r, \theta) @f$ coordinate.
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
        // The function describing how the derivative of the evolve function is calculated.
        auto dy = [&](DTensor<PseudoCartBasis>& updated_advection_field, CoordRTheta const& foot) {
            DVector<R, Theta> updated_advection_field_adv_space;
            ddcHelper::get<R>(updated_advection_field_adv_space) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<R>(m_advection_field_coefs)[idx_batch]));
            ddcHelper::get<Theta>(updated_advection_field_adv_space) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<Theta>(m_advection_field_coefs)[idx_batch]));
            // Ensure coord is inside the domain as splines can't extrapolate
            // derivates (clamping)
            CoordRTheta advection_location_for_mapping(
                    Kokkos::min(ddc::select<R>(foot), ddc::discrete_space<BSplinesR>().rmax()),
                    ddc::select<Theta>(foot));
            updated_advection_field = to_vector_space<PseudoCartBasis>(
                    m_logical_to_pseudo_physical,
                    advection_location_for_mapping,
                    updated_advection_field_adv_space);
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DTensor<PseudoCartBasis> const& advection_field,
                                   double dt) {
            CoordX_pcY_pc const coord_xy = m_logical_to_pseudo_physical(foot_rtheta);
            CoordX_pcY_pc const foot_xy = coord_xy - dt * advection_field;

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

/**
 * @brief Owning manager for @ref ElementwiseLogicalAdvPseudoPhysFootFinder.
 *
 * Holds the spline coefficient field and the preallocated time stepper. Call
 * @c operator()(dt) to obtain a GPU-copyable
 * @ref ElementwiseLogicalAdvPseudoPhysFootFinder configured for time step @p dt.
 */
template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefFieldMem,
        class TimeStepperBuilder,
        class LogicalToPseudoPhysicalMapping>
class ElementwiseLogicalAdvPseudoPhysFootFinderMem
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;
    using CoordRTheta = Coord<R, Theta>;
    using PseudoPhysicalToLogicalMapping = inverse_mapping_t<LogicalToPseudoPhysicalMapping>;
    using CoordX_pcY_pc = typename PseudoPhysicalToLogicalMapping::CoordArg;
    using PseudoCartBasis = ddc::to_type_seq_t<CoordX_pcY_pc>;

    using TimeStepper = typename TimeStepperBuilder::
            template time_stepper_t<CoordRTheta, DTensor<PseudoCartBasis>>;

public:
    /// The non-owning operator that can be used on GPU.
    using GPUCompat = ElementwiseLogicalAdvPseudoPhysFootFinder<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            typename AdvecCoefFieldMem::view_type,
            TimeStepper,
            LogicalToPseudoPhysicalMapping>;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;
    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;
    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    TimeStepper m_time_stepper;
    AdvecCoefFieldMem m_advection_field_coefs_alloc;
    CoordX_pcY_pc m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;

public:
    /**
     * @brief Construct an ElementwiseLogicalAdvPseudoPhysFootFinderMem.
     *
     * @param[in] evaluator_advection_field
     *      The evaluator for the spline representation of the advection field.
     * @param[in] logical_to_physical
     *      The mapping from the logical domain to the physical domain, used to
     *      derive the pseudo-physical mappings.
     * @param[in] time_stepper_builder
     *      The factory used to preallocate the time integration method.
     * @param[in] advection_field_coefs
     *      The spline coefficients of the advection field. Ownership is transferred in.
     * @param[in] coord_centre
     *      The coordinate of the polar centre in the pseudo-physical domain,
     *      used to handle the degenerate point at @f$ r = 0 @f$.
     * @param[in] idx_range_theta
     *      The poloidal index range, used to wrap the angular coordinate into
     *      the periodic domain after each time step.
     */
    template <
            class LogicalToPhysicalMapping,
            std::enable_if_t<
                    ddc::type_seq_contains_v<
                            PseudoCartBasis,
                            ddc::to_type_seq_t<typename Mapping::CoordResult>>,
                    bool> = true>
    ElementwiseLogicalAdvPseudoPhysFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            CoordX_pcY_pc coord_centre,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_logical_to_pseudo_physical(logical_to_physical)
        , m_pseudo_physical_to_logical(logical_to_physical.get_inverse_mapping())
        , m_time_stepper(time_stepper_builder.template preallocate<TimeStepper>())
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_coord_centre(coord_centre)
        , m_idx_range_theta(idx_range_theta)
    {
    }
    /**
     * @brief Construct an ElementwiseLogicalAdvPseudoPhysFootFinderMem when the
     *        logical-to-physical mapping does not share the pseudo-physical Cartesian tags.
     *
     * In this case @p logical_to_physical is ignored and the pseudo-physical mapping is
     * constructed directly from @p coord_centre using a default circular mapping.
     *
     * @param[in] evaluator_advection_field
     *      The evaluator for the spline representation of the advection field.
     * @param[in] logical_to_physical
     *      Unused; accepted for interface uniformity with the primary constructor.
     * @param[in] time_stepper_builder
     *      The factory used to preallocate the time integration method.
     * @param[in] advection_field_coefs
     *      The spline coefficients of the advection field. Ownership is transferred in.
     * @param[in] coord_centre
     *      The coordinate of the polar centre in the pseudo-physical domain,
     *      used to initialise the pseudo-physical mappings and handle the degenerate
     *      point at @f$ r = 0 @f$.
     * @param[in] idx_range_theta
     *      The poloidal index range, used to wrap the angular coordinate into
     *      the periodic domain after each time step.
     */
    template <
            class LogicalToPhysicalMapping,
            std::enable_if_t<
                    !(ddc::type_seq_contains_v<
                            PseudoCartBasis,
                            ddc::to_type_seq_t<typename Mapping::CoordResult>>),
                    bool> = false>
    ElementwiseLogicalAdvPseudoPhysFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            [[maybe_unused]] LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            CoordX_pcY_pc coord_centre,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_logical_to_pseudo_physical(coord_centre)
        , m_pseudo_physical_to_logical(coord_centre)
        , m_time_stepper(time_stepper_builder.template preallocate<TimeStepper>())
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_coord_centre(coord_centre)
        , m_idx_range_theta(idx_range_theta)
    {
    }

    /**
     * @brief Create an ElementwiseLogicalAdvPseudoPhysFootFinder for the given time step.
     *
     * Returns a non-owning @ref ElementwiseLogicalAdvPseudoPhysFootFinder that holds views of
     * the stored spline coefficients and is configured for time step @f$ dt @f$.
     * The returned object can be copied to and called from the device.
     *
     * @param[in] dt
     *      The time step for the characteristic equation.
     *
     * @return A view-based ElementwiseLogicalAdvPseudoPhysFootFinder for the given time step.
     */
    GPUCompat operator()(double dt)
    {
        return GPUCompat(
                m_evaluator_advection_field,
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
 * @brief Selects @ref ElementwiseLogicalAdvPseudoPhysFootFinderMem (circular pseudo-physical
 *        mapping) for logical advection with foot-finding in pseudo-physical space.
 */
template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping>
struct ElementwiseChoice<
        FootFindingSpace::PSEUDO_PHYSICAL,
        AdvectionFieldSpace::LOGICAL,
        GridR,
        GridTheta,
        IdxRangeOperator,
        RThetaAdvectionEvaluator,
        AdvecCoefField,
        TimeStepperBuilder,
        LogicalToPhysicalMapping>
{
    /// The selected elementwise operator type.
    using type = ElementwiseLogicalAdvPseudoPhysFootFinderMem<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            AdvecCoefField,
            TimeStepperBuilder,
            CircularToCartesian<
                    typename GridR::continuous_dimension_type,
                    typename GridTheta::continuous_dimension_type,
                    X_pC,
                    Y_pC>>;
};

/**
 * @brief Selects @ref ElementwiseLogicalAdvPseudoPhysFootFinderMem (using the physical
 *        mapping directly) for logical advection with foot-finding in physical space.
 */
template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping>
struct ElementwiseChoice<
        FootFindingSpace::PHYSICAL,
        AdvectionFieldSpace::LOGICAL,
        GridR,
        GridTheta,
        IdxRangeOperator,
        RThetaAdvectionEvaluator,
        AdvecCoefField,
        TimeStepperBuilder,
        LogicalToPhysicalMapping>
{
    static_assert(is_analytical_mapping_v<LogicalToPhysicalMapping>);
    /// The selected elementwise operator type.
    using type = ElementwiseLogicalAdvPseudoPhysFootFinderMem<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            AdvecCoefField,
            TimeStepperBuilder,
            LogicalToPhysicalMapping>;
};

} // namespace polar_foot_finder_details
