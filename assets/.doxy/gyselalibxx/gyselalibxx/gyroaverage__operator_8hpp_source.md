

# File gyroaverage\_operator.hpp

[**File List**](files.md) **>** [**gyroaverage**](dir_0488de83d5a156211532c59817863a53.md) **>** [**gyroaverage\_operator.hpp**](gyroaverage__operator_8hpp.md)

[Go to the documentation of this file](gyroaverage__operator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "i_interpolation.hpp"

template <
        concepts::Interpolation RThetaInterpolator,
        class IdxRangeRminorThetaBatch,
        class ToLogicalCoordTransform>
class GyroAverageOperator
{
    using RThetaBuilder = typename RThetaInterpolator::BuilderType;
    using RThetaEvaluator = typename RThetaInterpolator::EvaluatorType;

    using ExecutionSpace = typename RThetaBuilder::exec_space;

    using Rminor = typename CoordWithOPoint<
            typename ToLogicalCoordTransform::CoordResult>::curvilinear_tag_r;
    using Theta = typename CoordWithOPoint<
            typename ToLogicalCoordTransform::CoordResult>::curvilinear_tag_theta;

    using IdxRangeRminorTheta
            = InterpolationBuilderTraits<RThetaBuilder>::interpolation_idx_range_type;

    using GridRminor = find_grid_t<Rminor, ddc::to_type_seq_t<IdxRangeRminorTheta>>;
    using GridTheta = find_grid_t<Theta, ddc::to_type_seq_t<IdxRangeRminorTheta>>;

    struct R_gyro_cov;
    struct Theta_gyro_cov;

    struct R_gyro
    {
        static constexpr bool PERIODIC = false;
        // The corresponding type in the dual space.
        using Dual = R_gyro_cov;
    };

    struct Theta_gyro
    {
        static constexpr bool PERIODIC = true;
        // The corresponding type in the dual space.
        using Dual = Theta_gyro_cov;
    };

    using IdxRangeRminor = IdxRange<GridRminor>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeRminorThetaBatch, GridRminor, GridTheta>;
    using IdxRangeBSRminorTheta = InterpolationBuilderTraits<RThetaBuilder>::coeff_idx_range_type;

    using IdxRminor = Idx<GridRminor>;
    using IdxTheta = Idx<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxRminorTheta = Idx<GridRminor, GridTheta>;

    using DFieldMemRminorTheta = DFieldMem<IdxRangeRminorTheta>;
    using DFieldMemRminorThetaBatch = DFieldMem<IdxRangeRminorThetaBatch>;
    using DFieldMemBSRminorTheta = DFieldMem<IdxRangeBSRminorTheta>;
    using DFieldRminorTheta = DField<IdxRangeRminorTheta>;
    using DFieldRminorThetaBatch = DField<IdxRangeRminorThetaBatch>;
    using DFieldBSRminorTheta = DField<IdxRangeBSRminorTheta>;
    using DConstFieldRminorTheta = DConstField<IdxRangeRminorTheta>;
    using DConstFieldRminorThetaBatch = DConstField<IdxRangeRminorThetaBatch>;

    using CoordRminorTheta = Coord<Rminor, Theta>;
    using CoordR_gyroTheta_gyro = Coord<R_gyro, Theta_gyro>;
    using RZ_1 = ddc::
            type_seq_element_t<0, ddc::to_type_seq_t<typename ToLogicalCoordTransform::CoordArg>>;
    using RZ_2 = ddc::
            type_seq_element_t<1, ddc::to_type_seq_t<typename ToLogicalCoordTransform::CoordArg>>;
    using CoordRZ = typename ToLogicalCoordTransform::CoordArg;

    // FIXME
    // Need to add a static assert to check evaluator is addmissible to builder
    // using ddc::is_evaluator_admissible
    // This will be available from DDC 0.9.0
    static_assert(
            is_mapping_v<ToLogicalCoordTransform>,
            "CoordinateTransformFunction must be a mapping");
    static_assert(std::is_same_v<typename ToLogicalCoordTransform::CoordArg, CoordRZ>);
    static_assert(std::is_same_v<typename ToLogicalCoordTransform::CoordResult, CoordRminorTheta>);
    static_assert(is_accessible_v<ExecutionSpace, ToLogicalCoordTransform>);

    DConstFieldRminorTheta m_rho_L;

    RThetaBuilder const& m_builder;

    RThetaEvaluator const& m_evaluator;

    ToLogicalCoordTransform m_coordinate_transform;

    std::size_t const m_nb_gyro_points;

public:
    explicit GyroAverageOperator(
            DConstFieldRminorTheta const& rho_L,
            RThetaInterpolator const& interpolator,
            ToLogicalCoordTransform coordinate_transform,
            std::size_t const nb_gyro_points = 8)
        : m_rho_L(rho_L)
        , m_builder(interpolator.get_builder())
        , m_evaluator(interpolator.get_evaluator())
        , m_coordinate_transform(coordinate_transform)
        , m_nb_gyro_points(nb_gyro_points)
    {
    }

    void operator()(DFieldRminorThetaBatch const& A_bar, DConstFieldRminorThetaBatch const& A) const
    {
        IdxRangeRminorThetaBatch const rthetabatch_idx_range = get_idx_range(A);
        IdxRangeTheta const theta_idx_range(rthetabatch_idx_range);
        IdxRangeBatch const batch_idx_range(rthetabatch_idx_range);
        IdxRangeRminorTheta const rtheta_idx_range(rthetabatch_idx_range);

        // Instantiate chunk of interpolation coefs to receive output of builder (r, theta)
        DFieldMemBSRminorTheta coef_alloc(
                "coef (GyroAvrageOperator::operator())",
                batched_basis_idx_range(m_builder, rtheta_idx_range));
        DFieldBSRminorTheta const coef = get_field(coef_alloc);
        DConstFieldRminorTheta const rho_L = get_const_field(m_rho_L);

        using SubConstDFieldRminorTheta = DConstField<
                IdxRangeRminorTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;
        using SubDFieldRminorTheta = DField<
                IdxRangeRminorTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;

        ddc::host_for_each(batch_idx_range, [&](IdxBatch const ib) {
            SubConstDFieldRminorTheta const sub_A = A[ib];
            SubDFieldRminorTheta sub_A_bar = A_bar[ib];
            DFieldMemRminorTheta
                    sub_A_alloc("sub_A (GyroAvrageOperator::operator())", rtheta_idx_range);
            ddc::parallel_deepcopy(sub_A_alloc, sub_A);
            m_builder(coef, get_const_field(sub_A_alloc));
            RThetaEvaluator evaluator = m_evaluator;

            ToLogicalCoordTransform coordinate_transform = m_coordinate_transform;
            std::size_t nb_gyro_points = m_nb_gyro_points;
            // Loop over r, theta
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    ExecutionSpace(),
                    rtheta_idx_range,
                    KOKKOS_LAMBDA(IdxRminorTheta const irtheta) {
                        IdxRminor const ir(irtheta);
                        IdxTheta const itheta(irtheta);

                        // Average over gyro points
                        double sum_over_gyro_points = 0.0;
                        for (std::size_t igyro = 0; igyro < nb_gyro_points; igyro++) {
                            // Compute the particle position in (R, Z) coordinate
                            double const alpha = M_PI * 2.0 / static_cast<double>(nb_gyro_points)
                                                 * static_cast<double>(igyro);
                            inverse_mapping_t<ToLogicalCoordTransform> inv_coordinate_transform
                                    = coordinate_transform.get_inverse_mapping();
                            CoordRZ gyrocentre = inv_coordinate_transform(ddc::coordinate(irtheta));
                            CircularToCartesian<R_gyro, Theta_gyro, RZ_1, RZ_2> circ_to_cart(
                                    gyrocentre);
                            CoordRZ particle_position = circ_to_cart(
                                    CoordR_gyroTheta_gyro {rho_L(ir, itheta), alpha});

                            // Convert from (R, Z) into (r, theta) coordinate
                            CoordRminorTheta p = coordinate_transform(particle_position);

                            // Interpolation in (r, theta) coordinate
                            sum_over_gyro_points += evaluator(p, get_const_field(coef));
                        }
                        sub_A_bar(ir, itheta)
                                = sum_over_gyro_points / static_cast<double>(nb_gyro_points);
                    });

            // Apply periodic boundary condition in theta direction
            ddc::parallel_deepcopy(
                    sub_A_bar[theta_idx_range.back()],
                    sub_A_bar[theta_idx_range.front()]);
        });
    }
};
```


