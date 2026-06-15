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
        class TimeStepper,
        class LogicalToPhysicalMapping>
class ElementwisePhysicalAdvPhysicalFootFinder
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

    using PhysicalToLogicalMapping = inverse_mapping_t<LogicalToPhysicalMapping>;

    using VectorIndexSetAdvectionDims
            = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;
    /// The first dimension of the advection field.
    using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
    /// The second dimension of the advection field.
    using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;

    TimeStepper m_time_stepper;

    LogicalToPhysicalMapping m_logical_to_physical;
    PhysicalToLogicalMapping m_physical_to_logical;

    AdvecCoefField m_advection_field_coefs;
    IdxRange<GridTheta> m_idx_range_theta;
    Coord<AdvDim1, AdvDim2> m_coord_centre;
    double m_dt;

public:
    ElementwisePhysicalAdvPhysicalFootFinder(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            TimeStepper const& time_stepper,
            LogicalToPhysicalMapping logical_to_physical,
            AdvecCoefField const& advection_field_coefs,
            IdxRange<GridTheta> idx_range_theta,
            Coord<AdvDim1, AdvDim2> coord_centre,
            double dt)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_time_stepper(time_stepper)
        , m_logical_to_physical(logical_to_physical)
        , m_physical_to_logical(logical_to_physical.get_inverse_mapping())
        , m_advection_field_coefs(advection_field_coefs)
        , m_idx_range_theta(idx_range_theta)
        , m_coord_centre(coord_centre)
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
        auto dy = [&](DVector<AdvDim1, AdvDim2>& updated_advection_field, CoordRTheta const& foot) {
            ddcHelper::get<AdvDim1>(updated_advection_field) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<AdvDim1>(m_advection_field_coefs)[idx_batch]));
            ddcHelper::get<AdvDim2>(updated_advection_field) = m_evaluator_advection_field(
                    foot,
                    get_const_field(ddcHelper::get<AdvDim2>(m_advection_field_coefs)[idx_batch]));
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DVector<AdvDim1, AdvDim2> const& advection_field,
                                   double dt) {
            Coord<AdvDim1, AdvDim2> const coord_xy = m_logical_to_physical(foot_rtheta);
            Coord<AdvDim1, AdvDim2> const foot_xy = coord_xy - dt * advection_field;

            if (norm_inf(foot_xy - m_coord_centre) < 1e-15) {
                foot_rtheta = CoordRTheta(0, 0);
            } else {
                foot_rtheta = m_physical_to_logical(foot_xy);
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
        class LogicalToPhysicalMapping>
class ElementwisePhysicalAdvPhysicalFootFinderMem
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;
    using CoordRTheta = Coord<R, Theta>;

    using VectorIndexSetAdvectionDims
            = ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;
    /// The first dimension of the advection field.
    using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
    /// The second dimension of the advection field.
    using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

    using TimeStepper = typename TimeStepperBuilder::
            template time_stepper_t<CoordRTheta, DVector<AdvDim1, AdvDim2>>;

public:
    /// The non-owning operator that can be used on GPU.
    using GPUCompat = ElementwisePhysicalAdvPhysicalFootFinder<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            typename AdvecCoefFieldMem::view_type,
            TimeStepper,
            LogicalToPhysicalMapping>;

private:
    RThetaAdvectionEvaluator m_evaluator_advection_field;
    TimeStepper m_time_stepper;
    LogicalToPhysicalMapping m_logical_to_physical;
    AdvecCoefFieldMem m_advection_field_coefs_alloc;
    IdxRange<GridTheta> m_idx_range_theta;
    Coord<AdvDim1, AdvDim2> m_coord_centre;

public:
    /**
     * @brief Construct an ElementwisePhysicalAdvPhysicalFootFinderMem.
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
    ElementwisePhysicalAdvPhysicalFootFinderMem(
            RThetaAdvectionEvaluator const& evaluator_advection_field,
            LogicalToPhysicalMapping const& logical_to_physical,
            TimeStepperBuilder const& time_stepper_builder,
            AdvecCoefFieldMem&& advection_field_coefs,
            Coord<AdvDim1, AdvDim2> coord_centre,
            IdxRange<GridTheta> idx_range_theta)
        : m_evaluator_advection_field(evaluator_advection_field)
        , m_time_stepper(time_stepper_builder.template preallocate<TimeStepper>())
        , m_logical_to_physical(logical_to_physical)
        , m_advection_field_coefs_alloc(std::move(advection_field_coefs))
        , m_idx_range_theta(idx_range_theta)
        , m_coord_centre(coord_centre)
    {
    }

    /**
     * @brief Create an ElementwisePhysicalAdvPhysicalFootFinder for the given time step.
     *
     * Returns a non-owning @ref ElementwisePhysicalAdvPhysicalFootFinder that holds views of
     * the stored spline coefficients and is configured for time step @f$ dt @f$.
     * The returned object can be copied to and called from the device.
     *
     * @param[in] dt
     *      The time step for the characteristic equation.
     *
     * @return A view-based ElementwisePhysicalAdvPhysicalFootFinder for the given time step.
     */
    GPUCompat operator()(double dt)
    {
        return GPUCompat(
                m_evaluator_advection_field,
                m_time_stepper,
                m_logical_to_physical,
                get_const_field(m_advection_field_coefs_alloc),
                m_idx_range_theta,
                m_coord_centre,
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
        FootFindingSpace::PHYSICAL,
        AdvectionFieldSpace::PHYSICAL,
        GridR,
        GridTheta,
        IdxRangeOperator,
        RThetaAdvectionEvaluator,
        AdvecCoefField,
        TimeStepperBuilder,
        LogicalToPhysicalMapping>
{
    using type = ElementwisePhysicalAdvPhysicalFootFinderMem<
            GridR,
            GridTheta,
            IdxRangeOperator,
            RThetaAdvectionEvaluator,
            AdvecCoefField,
            TimeStepperBuilder,
            LogicalToPhysicalMapping>;
};

} // namespace polar_foot_finder_details
