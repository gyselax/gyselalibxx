// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"

namespace detail {

/**
 * @brief Operator to compute the gyroaverage of a field in (r, theta) coordinates.
 *
 * This class performs the gyroaveraging operation on a batched field defined on a polar grid (r, theta).
 * The gyroaverage is computed by integrating the field over a set of points along a circle (the Larmor orbit)
 * centred at each grid point, with radius given by the local Larmor radius field (rho_L).
 *
 * The class uses 2D B-spline interpolation to evaluate the field at off-grid points along the orbit.
 * The operation is performed in parallel over the (r, theta) grid.
 *
 * @tparam SplineRThetaBuilder      The type of the spline builder for the rtheta interpolation
 * @tparam SplineRThetaEvaluator    The type of the spline evaluator for the rtheta interpolation
 * @tparam IdxRangeRminorThetaBatch The index range over R, Theta and Batch directions.
 * @tparam CoordinateTransformFunction The Coordinate transform Functor.
 */
template <
        class SplineRThetaBuilder,
        class SplineRThetaEvaluator,
        class IdxRangeRminorThetaBatch,
        class CoordinateTransformFunction>
class GyroAverageOperator
{
    // FIXME
    // Need to add a static assert to check evaluator is addmissible to builder
    // using ddc::is_evaluator_admissible
    // This will be available from DDC 0.9.0
    using ExecutionSpace = typename SplineRThetaBuilder::exec_space;

    using GridRminor = typename SplineRThetaBuilder::interpolation_discrete_dimension_type1;
    using GridTheta = typename SplineRThetaBuilder::interpolation_discrete_dimension_type2;

    using Rminor = typename GridRminor::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    using BSplinesRminor = typename SplineRThetaBuilder::bsplines_type1;
    using BSplinesTheta = typename SplineRThetaBuilder::bsplines_type2;

public:
    using IdxRangeRminor = IdxRange<GridRminor>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeRminorTheta = IdxRange<GridRminor, GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeRminorThetaBatch, GridRminor, GridTheta>;
    using IdxRangeBSRminorTheta = IdxRange<BSplinesRminor, BSplinesTheta>;

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
    using DConstFieldRminorTheta = ConstField<double, IdxRangeRminorTheta>;
    using DConstFieldRminorThetaBatch = ConstField<double, IdxRangeRminorThetaBatch>;

    using CoordRminorTheta = Coord<Rminor, Theta>;
    using CoordRZ = CoordXY_pC;

private:
    /**
     * @brief Field of local Larmor radii (rho_L) on the (r, theta) grid.
     */
    DConstFieldRminorTheta m_rho_L;

    /**
     * @brief The spline builder for the rtheta interpolation
     */
    SplineRThetaBuilder const& m_spline_builder;

    /**
     * @brief The spline evaluator for the rtheta interpolation
     */
    SplineRThetaEvaluator const& m_spline_evaluator;

    /**
     * @brief coordinate_transform Function to convert (R, Z) to (r, theta).
     */
    CoordinateTransformFunction m_coordinate_transform;

    /**
     * @brief Number of points to use in the gyroaverage integration (default: 8)
     */
    std::size_t const m_nb_gyro_points;

public:
    /**
     * @brief Constructor.
     * @param[in] rho_L Field of Larmor radii on the (r, theta) grid.
     * @param[in] spline_builder The spline builder for the rtheta interpolation
     * @param[in] spline_evaluator The spline evaluator for the rtheta interpolation
     * @param[in] coordinate_transform Function to convert (R, Z) to (r, theta).
     * @param[in] nb_gyro_points Number of points to use in the gyroaverage integration (default: 8).
     */
    explicit GyroAverageOperator(
            DConstFieldRminorTheta const& rho_L,
            SplineRThetaBuilder const& spline_builder,
            SplineRThetaEvaluator const& spline_evaluator,
            CoordinateTransformFunction coordinate_transform,
            std::size_t const nb_gyro_points = 8)
        : m_rho_L(rho_L)
        , m_spline_builder(spline_builder)
        , m_spline_evaluator(spline_evaluator)
        , m_coordinate_transform(coordinate_transform)
        , m_nb_gyro_points(nb_gyro_points)
    {
        static_assert(
                is_mapping_v<CoordinateTransformFunction>,
                "CoordinateTransformFunction must be a mapping");
        static_assert(std::is_same_v<typename CoordinateTransformFunction::CoordArg, CoordRZ>);
        static_assert(std::is_same_v<
                      typename CoordinateTransformFunction::CoordResult,
                      CoordRminorTheta>);
        static_assert(is_accessible_v<ExecutionSpace, CoordinateTransformFunction>);
    }

    /**
     * @brief Applies the gyroaverage operator to a batched field.
     *
     * For each batch, and for each (r, theta) grid point, computes the gyroaverage by
     * integrating the field along a circle of radius rho_L centred at (r, theta).
     * The field is interpolated at off-grid points using 2D B-splines.
     *
     * @tparam CoordinateTransformFunction
     *         Callable type that transforms (R, Z) coordinates to (r, theta) coordinates.
     * @param[in] A Input field to be gyroaveraged (batched).
     * @param[out] A_bar Output field to store the gyroaveraged result (batched).
     * @param[in] coordinate_transform Function to convert (R, Z) to (r, theta).
     */
    void operator()(DConstFieldRminorThetaBatch const& A, DFieldRminorThetaBatch const& A_bar) const
    {
        IdxRangeRminorThetaBatch const rthetabatch_idx_range = get_idx_range(A);
        IdxRangeTheta const theta_idx_range(rthetabatch_idx_range);
        IdxRangeBatch const batch_idx_range(rthetabatch_idx_range);
        IdxRangeRminorTheta const rtheta_idx_range(rthetabatch_idx_range);

        // Instantiate chunk of spline coefs to receive output of spline_builder (r, theta)
        DFieldMemBSRminorTheta coef_alloc(get_spline_idx_range(m_spline_builder));
        DFieldBSRminorTheta const coef = get_field(coef_alloc);
        DConstFieldRminorTheta const rho_L = get_const_field(m_rho_L);

        using SubConstDFieldRminorTheta = ConstField<
                double,
                IdxRangeRminorTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;
        using SubDFieldRminorTheta = DField<
                IdxRangeRminorTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;

        ddc::for_each(batch_idx_range, [&](IdxBatch const ib) {
            // FIXME
            // The input of the spline builder must be LayoutRight
            // We allocate a buffer in LayoutRight whereto the slice is copied
            SubConstDFieldRminorTheta const sub_A = A[ib];
            SubDFieldRminorTheta sub_A_bar = A_bar[ib];
            DFieldMemRminorTheta sub_A_alloc(rtheta_idx_range);
            ddc::parallel_deepcopy(sub_A_alloc, sub_A);
            m_spline_builder(coef, get_const_field(sub_A_alloc));
            SplineRThetaEvaluator spline_evaluator = m_spline_evaluator;

            CoordinateTransformFunction coordinate_transform = m_coordinate_transform;
            std::size_t nb_gyro_points = m_nb_gyro_points;
            // Loop over r, theta
            ddc::parallel_for_each(
                    ExecutionSpace(),
                    rtheta_idx_range,
                    KOKKOS_LAMBDA(IdxRminorTheta const irtheta) {
                        IdxRminor const ir(irtheta);
                        IdxTheta const itheta(irtheta);
                        double const r = ddc::coordinate(ir);
                        double const theta = ddc::coordinate(itheta);

                        // Average over gyro points
                        double sum_over_gyro_points = 0.0;
                        for (std::size_t igyro = 0; igyro < nb_gyro_points; igyro++) {
                            // Compute the particle position in (R, Z) coordinate
                            double const alpha = M_PI * 2.0 / static_cast<double>(nb_gyro_points)
                                                 * static_cast<double>(igyro);
                            double const R_p = r * Kokkos::cos(theta)
                                               + rho_L(ir, itheta) * Kokkos::cos(alpha);
                            double const Z_p = r * Kokkos::sin(theta)
                                               + rho_L(ir, itheta) * Kokkos::sin(alpha);

                            // Convert from (R, Z) into (r, theta) coordinate
                            CoordRminorTheta p = coordinate_transform(CoordRZ {R_p, Z_p});

                            // Spline interpolation in (r, theta) coordinate
                            sum_over_gyro_points += spline_evaluator(p, get_const_field(coef));
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

} // namespace detail
