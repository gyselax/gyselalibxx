

# File spline\_polar\_foot\_finder.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**spline\_polar\_foot\_finder.hpp**](spline__polar__foot__finder_8hpp.md)

[Go to the documentation of this file](spline__polar__foot__finder_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "ipolar_foot_finder.hpp"
#include "mapping_tools.hpp"
#include "vector_index_tools.hpp"
#include "vector_mapper.hpp"

template <
        class TimeStepper,
        class LogicalToPhysicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class SplineRThetaBuilderAdvection,
        class SplineRThetaEvaluatorAdvection>
class SplinePolarFootFinder
    : public IPolarFootFinder<
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
              typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
              ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
              typename TimeStepper::IdxRange,
              typename SplineRThetaBuilderAdvection::memory_space>
{
    static_assert(is_mapping_v<LogicalToPhysicalMapping>);
    static_assert(is_mapping_v<LogicalToPseudoPhysicalMapping>);
    static_assert(is_analytical_mapping_v<LogicalToPseudoPhysicalMapping>);
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
                  ddc::to_type_seq_t<typename TimeStepper::IdxRange>>);
    static_assert(ddc::in_tags_v<
                  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
                  ddc::to_type_seq_t<typename TimeStepper::IdxRange>>);

    using base_type = IPolarFootFinder<
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1,
            typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2,
            ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>,
            typename TimeStepper::IdxRange,
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

private:
    using ExecSpace = typename SplineRThetaBuilderAdvection::exec_space;
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
    using CoordXY_adv = typename LogicalToPseudoPhysicalMapping::CoordResult;

    using PseudoCartesianBasis = ddc::to_type_seq_t<CoordXY_adv>;

    using CoordRTheta = Coord<R, Theta>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
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

    TimeStepper const& m_time_stepper;

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
            TimeStepper const& time_stepper,
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            LogicalToPseudoPhysicalMapping const& logical_to_pseudo_physical_mapping,
            SplineRThetaBuilderAdvection const& builder_advection_field,
            SplineRThetaEvaluatorAdvection const& evaluator_advection_field,
            double epsilon = 1e-12)
        : m_time_stepper(time_stepper)
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
        VectorSplineCoeffsMem advection_field_in_adv_domain_coefs(
                get_spline_idx_range(m_builder_advection_field));

        // Compute the advection field in the advection domain.
        auto advection_field_in_adv_domain = create_geometry_mirror_view(
                ExecSpace(),
                advection_field,
                m_pseudo_physical_to_physical);

        // Get the coefficients of the advection field in the advection domain.
        m_builder_advection_field(
                ddcHelper::get<X_adv>(advection_field_in_adv_domain_coefs),
                ddcHelper::get<X_adv>(get_const_field(advection_field_in_adv_domain)));
        m_builder_advection_field(
                ddcHelper::get<Y_adv>(advection_field_in_adv_domain_coefs),
                ddcHelper::get<Y_adv>(get_const_field(advection_field_in_adv_domain)));


        // The function describing how the derivative of the evolve function is calculated.
        std::function<void(DVectorFieldAdvection, CConstFieldFeet)> dy
                = [&](DVectorFieldAdvection updated_advection_field, CConstFieldFeet feet) {
                      m_evaluator_advection_field(
                              get_field(ddcHelper::get<X_adv>(updated_advection_field)),
                              get_const_field(feet),
                              get_const_field(
                                      ddcHelper::get<X_adv>(advection_field_in_adv_domain_coefs)));
                      m_evaluator_advection_field(
                              get_field(ddcHelper::get<Y_adv>(updated_advection_field)),
                              get_const_field(feet),
                              get_const_field(
                                      ddcHelper::get<Y_adv>(advection_field_in_adv_domain_coefs)));
                  };

        CoordXY_adv coord_centre(m_logical_to_pseudo_physical(CoordRTheta(0, 0)));
        LogicalToPseudoPhysicalMapping logical_to_pseudo_physical_proxy
                = m_logical_to_pseudo_physical;
        PseudoPhysicalToLogicalMapping pseudo_physical_to_logical_proxy
                = m_pseudo_physical_to_logical;

        IdxRangeOperator idx_range = get_idx_range(feet);
        IdxRangeTheta idx_range_theta(idx_range);

        // The function describing how the value(s) are updated using the derivative.
        std::function<void(CFieldFeet, DVectorConstFieldAdvection, double)> update_function
                = [&](CFieldFeet feet, DVectorConstFieldAdvection advection_field, double dt) {
                      // Compute the characteristic feet at t^n:
                      ddc::parallel_for_each(
                              ExecSpace(),
                              get_idx_range(feet),
                              KOKKOS_LAMBDA(IdxOperator const idx) {
                                  CoordRTheta const coord_rtheta(feet(idx));
                                  CoordXY_adv const coord_xy
                                          = logical_to_pseudo_physical_proxy(coord_rtheta);

                                  CoordXY_adv const feet_xy = coord_xy - dt * advection_field(idx);

                                  if (norm_inf(feet_xy - coord_centre) < 1e-15) {
                                      feet(idx) = CoordRTheta(0, 0);
                                  } else {
                                      feet(idx) = pseudo_physical_to_logical_proxy(feet_xy);
                                      ddc::select<Theta>(feet(idx))
                                              = ddcHelper::restrict_to_idx_range(
                                                      ddc::select<Theta>(feet(idx)),
                                                      idx_range_theta);
                                  }
                              });

                      // Treatment to conserve the C0 property of the advected function:
                      unify_value_at_centre_pt(feet);
                      // Test if the values are the same at the centre point
                      is_unified(feet);
                  };


        // Solve the characteristic equation
        m_time_stepper.update(ExecSpace(), feet, dt, dy, update_function);

        is_unified(feet);
    }



    template <class T>
    void is_unified(Field<T, IdxRangeRTheta, memory_space> const& values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        IdxR r0_idx = r_idx_range.front();
        if (Kokkos::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            ddc::parallel_for_each(
                    ExecSpace(),
                    theta_idx_range,
                    KOKKOS_LAMBDA(const IdxTheta itheta) {
                        if (norm_inf(
                                    values(r0_idx, itheta)
                                    - values(r0_idx, theta_idx_range.front()))
                            > 1e-15) {
                            Kokkos::printf("WARNING ! -> Discontinuous at the centre point.");
                        }
                        KOKKOS_ASSERT(
                                values(r0_idx, itheta) == values(r0_idx, theta_idx_range.front()));
                    });
        }
    }


    template <class T>
    void unify_value_at_centre_pt(Field<T, IdxRangeRTheta, memory_space> values) const
    {
        IdxRangeR const r_idx_range = get_idx_range<GridR>(values);
        IdxRangeTheta const theta_idx_range = get_idx_range<GridTheta>(values);
        IdxR r0_idx = r_idx_range.front();
        if (std::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            ddc::parallel_for_each(
                    ExecSpace(),
                    theta_idx_range,
                    KOKKOS_LAMBDA(const IdxTheta itheta) {
                        values(r0_idx, itheta) = values(r0_idx, theta_idx_range.front());
                    });
        }
    }
};
```


