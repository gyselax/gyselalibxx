

# File spline\_polar\_foot\_finder.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**spline\_polar\_foot\_finder.hpp**](spline__polar__foot__finder_8hpp.md)

[Go to the documentation of this file](spline__polar__foot__finder_8hpp.md)


```C++
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

template <
        class IdxRangeBatched,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping,
        concepts::AnalyticalMapping LogicalToPseudoPhysicalMapping,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplinePolarFootFinder
    : public IPolarFootFinder<
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
              ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
              IdxRangeBatched,
              typename SplineRThetaBuilderAdvection::memory_space>
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

    using base_type = IPolarFootFinder<
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
            ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
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
    using ExecSpace = typename SplineRThetaBuilderAdvection::exec_space;

private:
    using MemSpace = typename ExecSpace::memory_space;
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
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
    using CFieldFeet = Field<CoordRTheta, IdxRangeOperator, memory_space>;

    using CConstFieldFeet = ConstField<CoordRTheta, IdxRangeOperator, memory_space>;

    using DVectorFieldAdvection
            = DVectorField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;

    using DVectorConstFieldAdvection
            = DVectorConstField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;

    using VectorSplineCoeffsMem
            = DVectorFieldMem<IdxRangeSplineBatched, PseudoCartesianBasis, memory_space>;

public:
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

    void operator()(
            CFieldFeet feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt) const final
    {
        static_assert(ddc::type_seq_size_v<VectorIndexSetAdvectionDims> == 2);
        using AdvDim1 = ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
        using AdvDim2 = ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;

        DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>
                advection_field_coefs_alloc(m_builder_advection_field.batched_spline_domain(
                        get_idx_range(advection_field)));
        DVectorField<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>
                advection_field_coefs(advection_field_coefs_alloc);

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<AdvDim1>(advection_field_coefs),
                ddcHelper::get<AdvDim1>(get_const_field(advection_field)));
        m_builder_advection_field(
                ddcHelper::get<AdvDim2>(advection_field_coefs),
                ddcHelper::get<AdvDim2>(get_const_field(advection_field)));

        CoordXY_adv coord_centre(m_logical_to_pseudo_physical(CoordRTheta(0, 0)));
        LogicalToPseudoPhysicalMapping logical_to_pseudo_physical_proxy
                = m_logical_to_pseudo_physical;
        PseudoPhysicalToLogicalMapping pseudo_physical_to_logical_proxy
                = m_pseudo_physical_to_logical;

        IdxRangeOperator idx_range = get_idx_range(feet);
        IdxRangeTheta idx_range_theta(idx_range);

        TimeStepper time_stepper = m_time_stepper_builder.template preallocate<TimeStepper>();

        SplineRThetaEvaluatorAdvection const& evaluator_advection_field_proxy
                = m_evaluator_advection_field;
        PseudoPhysicalToPhysicalMapping pseudo_physical_to_physical = m_pseudo_physical_to_physical;

        // Compute the characteristic feet at t^n:
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                get_idx_range(feet),
                KOKKOS_LAMBDA(IdxOperator const idx) {
                    IdxBatch idx_batch(idx);
                    IdxRTheta idx_rtheta(idx);
                    // The function describing how the derivative of the evolve function is calculated.
                    auto dy = [&](DVector<X_adv, Y_adv>& updated_advection_field,
                                  CoordRTheta const& foot) {
                        DVector<AdvDim1, AdvDim2> updated_advection_field_adv_space;
                        ddcHelper::get<AdvDim1>(updated_advection_field_adv_space)
                                = evaluator_advection_field_proxy(
                                        foot,
                                        get_const_field(ddcHelper::get<AdvDim1>(
                                                advection_field_coefs)[idx_batch]));
                        ddcHelper::get<AdvDim2>(updated_advection_field_adv_space)
                                = evaluator_advection_field_proxy(
                                        foot,
                                        get_const_field(ddcHelper::get<AdvDim2>(
                                                advection_field_coefs)[idx_batch]));
                        // Ensure coord is inside the domain as splines can't extrapolate
                        // derivates (clamping)
                        CoordRTheta advection_location_for_mapping(
                                Kokkos::
                                        min(ddc::select<R>(foot),
                                            ddc::discrete_space<BSplinesR>().rmax()),
                                ddc::select<Theta>(foot));
                        updated_advection_field = to_vector_space<VectorIndexSet<X_adv, Y_adv>>(
                                pseudo_physical_to_physical,
                                advection_location_for_mapping,
                                updated_advection_field_adv_space);
                    };

                    // The function describing how the value(s) are updated using the derivative.
                    auto update_function = [&](CoordRTheta& foot_rtheta,
                                               DVector<X_adv, Y_adv> const& advection_field,
                                               double dt) {
                        CoordXY_adv const coord_xy = logical_to_pseudo_physical_proxy(foot_rtheta);
                        CoordXY_adv const foot_xy = coord_xy - dt * advection_field;

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
```


