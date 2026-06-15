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
    using X_pc = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_pc = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;

    PseudoPhysicalToLogicalMapping m_pseudo_physical_to_logical;
    LogicalToPseudoPhysicalMapping m_logical_to_pseudo_physical;

    TimeStepper m_time_stepper;

    AdvecCoefField m_advection_field_coefs;
    Coord<X_pc, Y_pc> m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;
    double m_dt;

public:
    ElementwiseLogicalAdvPseudoPhysFootFinder(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            PseudoPhysicalToLogicalMapping const& pseudo_physical_to_logical,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical,
            TimeStepper const& time_stepper,
            AdvecCoefField const& advection_field_coefs,
            Coord<X_pc, Y_pc> coord_centre,
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

    KOKKOS_FUNCTION CoordRTheta operator()(IdxOperator const idx) const
    {
        IdxBatch idx_batch(idx);
        IdxRTheta idx_rtheta(idx);
        // The function describing how the derivative of the evolve function is calculated.
        auto dy = [&](DVector<X_pc, Y_pc>& updated_advection_field, CoordRTheta const& foot) {
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
            updated_advection_field = to_vector_space<VectorIndexSet<X_pc, Y_pc>>(
                    m_pseudo_physical_to_logical,
                    advection_location_for_mapping,
                    updated_advection_field_adv_space);
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DVector<X_pc, Y_pc> const& advection_field,
                                   double dt) {
            Coord<X_pc, Y_pc> const coord_xy = m_logical_to_pseudo_physical(foot_rtheta);
            Coord<X_pc, Y_pc> const foot_xy = coord_xy - dt * advection_field;

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
    using X_pc = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_pc = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;

    using TimeStepper =
            typename TimeStepperBuilder::template time_stepper_t<CoordRTheta, DVector<X_pc, Y_pc>>;

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
    Coord<X_pc, Y_pc> m_coord_centre;
    IdxRange<GridTheta> m_idx_range_theta;

public:
    /**
     * @brief Construct an ElementwiseLogicalAdvPseudoPhysFootFinderMem.
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
                    (std::is_same_v<X_pc, typename LogicalToPhysicalMapping::cartesian_tag_x>)&&(
                            std::is_same_v<
                                    Y_pc,
                                    typename LogicalToPhysicalMapping::cartesian_tag_y>),
                    bool> = true>
    ElementwiseLogicalAdvPseudoPhysFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            Coord<X_pc, Y_pc> coord_centre,
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
    template <
            class LogicalToPhysicalMapping,
            std::enable_if_t<
                    !((std::is_same_v<X_pc, typename LogicalToPhysicalMapping::cartesian_tag_x>)&&(
                            std::is_same_v<
                                    Y_pc,
                                    typename LogicalToPhysicalMapping::cartesian_tag_y>)),
                    bool> = false>
    ElementwiseLogicalAdvPseudoPhysFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            [[maybe_unused]] LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            Coord<X_pc, Y_pc> coord_centre,
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
