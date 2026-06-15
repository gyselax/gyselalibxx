// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "l_norm_tools.hpp"
#include "tensor.hpp"
#include "vector_field.hpp"

namespace polar_foot_finder_details {

template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepper>
class ElementwiseLogicalAdvLogicalFootFinder
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
    RThetaAdvectionEvaluator m_evaluator_advection_field;

    TimeStepper m_time_stepper;

    AdvecCoefField m_advection_field_coefs;
    IdxRange<GridTheta> m_idx_range_theta;
    double m_dt;

public:
    ElementwiseLogicalAdvLogicalFootFinder(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
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
        auto dy = [&](DVector<R, Theta>& updated_advection_field, CoordRTheta const& foot) {
            ddcHelper::get<R>(updated_advection_field) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<R>(m_advection_field_coefs)[idx_batch]));
            ddcHelper::get<Theta>(updated_advection_field) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<Theta>(m_advection_field_coefs)[idx_batch]));
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DVector<R, Theta> const& advection_field,
                                   double dt) {
            foot_rtheta -= dt * advection_field;

            // O-point reflection: if r goes negative, reflect through the origin.
            double foot_r = ddc::get<R>(foot_rtheta);
            double foot_theta = ddc::get<Theta>(foot_rtheta);
            int negative_reflexion = static_cast<int>(foot_r < 0);
            foot_rtheta = CoordRTheta(
                    (1 - 2 * negative_reflexion) * foot_r,
                    foot_theta + M_PI * negative_reflexion);
            // Wrap theta into the periodic domain.
            ddc::select<Theta>(foot_rtheta) = ddcHelper::
                    restrict_to_idx_range(ddc::select<Theta>(foot_rtheta), m_idx_range_theta);
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
        class TimeStepperBuilder>
class ElementwiseLogicalAdvLogicalFootFinderMem
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;
    using CoordRTheta = Coord<R, Theta>;

    using TimeStepper =
            typename TimeStepperBuilder::template time_stepper_t<CoordRTheta, DVector<R, Theta>>;


public:
    /// The non-owning operator that can be used on GPU.
    using GPUCompat = ElementwiseLogicalAdvLogicalFootFinder<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            typename AdvecCoefFieldMem::view_type,
            TimeStepper>;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;
    TimeStepper m_time_stepper;
    AdvecCoefFieldMem m_advection_field_coefs_alloc;
    IdxRange<GridTheta> m_idx_range_theta;

public:
    /**
     * @brief Construct an ElementwiseLogicalAdvLogicalFootFinderMem.
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
    template <class LogicalToPhysicalMapping, class X_pc, class Y_pc>
    ElementwiseLogicalAdvLogicalFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            [[maybe_unused]] LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            [[maybe_unused]] Coord<X_pc, Y_pc> coord_centre,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_time_stepper(time_stepper_builder.template preallocate<TimeStepper>())
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_idx_range_theta(idx_range_theta)
    {
    }

    /**
     * @brief Create an ElementwiseLogicalAdvLogicalFootFinder for the given time step.
     *
     * Returns a non-owning @ref ElementwiseLogicalAdvLogicalFootFinder that holds views of
     * the stored spline coefficients and is configured for time step @f$ dt @f$.
     * The returned object can be copied to and called from the device.
     *
     * @param[in] dt
     *      The time step for the characteristic equation.
     *
     * @return A view-based ElementwiseLogicalAdvLogicalFootFinder for the given time step.
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

template <
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping>
struct ElementwiseChoice<
        FootFindingSpace::LOGICAL,
        AdvectionFieldSpace::LOGICAL,
        GridR,
        GridTheta,
        IdxRangeOperator,
        RThetaAdvectionEvaluator,
        AdvecCoefField,
        TimeStepperBuilder,
        LogicalToPhysicalMapping>
{
    using type = ElementwiseLogicalAdvLogicalFootFinderMem<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            AdvecCoefField,
            TimeStepperBuilder>;
};

} // namespace polar_foot_finder_details
