

# File logical\_advection\_logical\_foot\_finder.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**logical\_advection\_logical\_foot\_finder.hpp**](logical__advection__logical__foot__finder_8hpp.md)

[Go to the documentation of this file](logical__advection__logical__foot__finder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "elementwise_choice.hpp"
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
    using R = typename GridR::continuous_dimension_type;
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

    KOKKOS_FUNCTION CoordRTheta operator()(IdxOperator const idx) const
    {
        IdxBatch idx_batch(idx);
        IdxRTheta idx_rtheta(idx);

        // The function describing how the derivative of the evolve function is calculated.
        auto dy = [&](DVector<R, Theta>& updated_advection_field, CoordRTheta const& foot) {
            int radial_sign;
            CoordRTheta foot_rtheta = rectify_coord_o_point(radial_sign, foot);

            ddcHelper::get<R>(updated_advection_field)
                    = m_evaluator_advection_field(
                              foot_rtheta,
                              get_const_field(
                                      ddcHelper::get<R>(m_advection_field_coefs)[idx_batch]))
                      * radial_sign; // Direction is reversed if r<0
            ddcHelper::get<Theta>(updated_advection_field) = m_evaluator_advection_field(
                    foot_rtheta,
                    get_const_field(ddcHelper::get<Theta>(m_advection_field_coefs)[idx_batch]));
        };

        // The function describing how the value(s) are updated using the derivative.
        auto update_function = [&](CoordRTheta& foot_rtheta,
                                   DVector<R, Theta> const& advection_field,
                                   double dt) {
            foot_rtheta -= dt * advection_field;
            // Wrap theta into the periodic domain.
            ddc::get<Theta>(foot_rtheta) = ddcHelper::
                    restrict_to_idx_range(ddc::select<Theta>(foot_rtheta), m_idx_range_theta);
        };

        CoordRTheta foot = ddc::coordinate(idx_rtheta);
        KOKKOS_ASSERT(ddc::select<R>(foot) > 1e-13);
        // Solve the characteristic equation
        m_time_stepper.update(foot, m_dt, dy, update_function);

        return rectify_coord_o_point(foot);
    }

private:
    KOKKOS_FUNCTION CoordRTheta rectify_coord_o_point(int& radial_sign, CoordRTheta foot) const
    {
        // O-point reflection: if r goes negative, reflect through the origin.
        double foot_r = ddc::get<R>(foot);
        double foot_theta = ddc::get<Theta>(foot);
        // 1 if r is negative, 0 otherwise
        int negative_reflection = static_cast<int>(foot_r < 0);
        // -1 if r is negative, 1 otherwise
        radial_sign = 1 - 2 * negative_reflection;
        // Get the equivalent foot with r>0
        CoordRTheta foot_rtheta(radial_sign * foot_r, foot_theta + M_PI * negative_reflection);
        // Wrap theta into the periodic domain.
        ddc::get<Theta>(foot_rtheta) = ddcHelper::
                restrict_to_idx_range(ddc::select<Theta>(foot_rtheta), m_idx_range_theta);
        return foot_rtheta;
    }

    KOKKOS_FUNCTION CoordRTheta rectify_coord_o_point(CoordRTheta foot) const
    {
        int radial_sign;
        return rectify_coord_o_point(radial_sign, foot);
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
```


