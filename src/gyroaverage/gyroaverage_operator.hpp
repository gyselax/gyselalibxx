// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"

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
 * @tparam ExecutionSpace      The Kokkos execution space.
 * @tparam GridR               The radial grid type.
 * @tparam GridTheta           The poloidal angle grid type.
 * @tparam GridBatch           The batch grid type (for batched operations).
 * @tparam BSplinesR           The B-spline basis for the radial direction.
 * @tparam BSplinesTheta       The B-spline basis for the theta direction.
 */
template <
        class ExecutionSpace,
        class GridR,
        class GridTheta,
        class GridBatch,
        class BSplinesR,
        class BSplinesTheta>
class GyroAverageOperator
{
    /**
     * @brief R grid.
     */
    using GridRType = GridR;
    /**
     * @brief Theta grid.
     */
    using GridThetaType = GridTheta;
    /**
     * @brief Batch grid.
     */
    using GridBatchType = GridBatch;

    using DimensionRType = typename GridRType::continuous_dimension_type;
    using DimensionThetaType = typename GridThetaType::continuous_dimension_type;

    using SplineRThetaBuilder = ddc::SplineBuilder2D<
            ExecutionSpace,
            typename ExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::BoundCond::GREVILLE, // boundary at r=0
            ddc::BoundCond::GREVILLE, // boundary at rmax
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK>;

    using SplineRThetaEvaluatorNullBound = ddc::SplineEvaluator2D<
            ExecutionSpace,
            typename ExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::NullExtrapolationRule, // boundary at r=0
            ddc::NullExtrapolationRule, // boundary at rmax
            ddc::PeriodicExtrapolationRule<DimensionThetaType>,
            ddc::PeriodicExtrapolationRule<DimensionThetaType>>;

public:
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeRThetaBatch = IdxRange<GridR, GridTheta, GridBatch>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeRThetaBatch, GridR, GridTheta>;
    using IdxRangeBSRTheta = IdxRange<BSplinesR, BSplinesTheta>;

    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxRTheta = Idx<GridR, GridTheta>;

    using DFieldMemRTheta = DFieldMem<IdxRangeRTheta>;
    using DFieldMemRThetaBatch = DFieldMem<IdxRangeRThetaBatch>;
    using DFieldMemBSRTheta = DFieldMem<IdxRangeBSRTheta>;
    using DFieldRTheta = DField<IdxRangeRTheta>;
    using DFieldRThetaBatch = DField<IdxRangeRThetaBatch>;
    using DFieldBSRTheta = DField<IdxRangeBSRTheta>;
    using DConstFieldRTheta = ConstField<double, IdxRangeRTheta>;
    using DConstFieldRThetaBatch = ConstField<double, IdxRangeRThetaBatch>;

    using CoordRTheta = ddc::Coordinate<DimensionRType, DimensionThetaType>;

private:
    /**
     * @brief Field of local Larmor radii (rho_L) on the (r, theta) grid.
     */
    DFieldRTheta m_rho_L;

public:
    /**
     * @brief Constructor.
     * @param[in] rho_L Field of Larmor radii on the (r, theta) grid.
     */
    explicit GyroAverageOperator(DFieldRTheta const& rho_L) : m_rho_L(rho_L) {}

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
     * @param[in] nb_gyro_points Number of points to use in the gyroaverage integration (default: 8).
     */
    template <class CoordinateTransformFunction>
    void operator()(
            DConstFieldRThetaBatch const& A,
            DFieldRThetaBatch const& A_bar,
            CoordinateTransformFunction coordinate_transform,
            std::size_t nb_gyro_points = 8) const
    {
        static_assert(
                std::is_invocable_v<CoordinateTransformFunction, double, double>,
                "CoordinateTransformFunction must be a functor on (double, double)");
        IdxRangeR const r_domain = get_idx_range<GridR>(A);
        IdxRangeTheta const theta_domain = get_idx_range<GridTheta>(A);
        IdxRangeBatch const batch_domain = get_idx_range<GridBatch>(A);
        IdxRangeRTheta const rtheta_mesh = get_idx_range<GridR, GridTheta>(A);

        ddc::NullExtrapolationRule r_extrapolation_rule;
        ddc::PeriodicExtrapolationRule<DimensionThetaType> theta_extrapolation_rule;

        SplineRThetaBuilder const spline_builder(rtheta_mesh);
        SplineRThetaEvaluatorNullBound const spline_evaluator(
                r_extrapolation_rule,
                r_extrapolation_rule,
                theta_extrapolation_rule,
                theta_extrapolation_rule);

        // Instantiate chunk of spline coefs to receive output of spline_builder (r, theta)
        DFieldMemBSRTheta coef_alloc(
                spline_builder.batched_spline_domain(rtheta_mesh),
                ddc::DeviceAllocator<double>());
        DFieldBSRTheta const coef = get_field(coef_alloc);
        DConstFieldRTheta const rho_L = get_const_field(m_rho_L);

        using SubConstDFieldRTheta = ConstField<
                double,
                IdxRangeRTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;
        using SubDFieldRTheta = DField<
                IdxRangeRTheta,
                typename ExecutionSpace::memory_space,
                Kokkos::layout_stride>;

        ddc::for_each(batch_domain, [&](IdxBatch const ib) {
            // FIXME
            // The input of the spline builder must be LayoutRight
            // We allocate a buffer in LayoutRight whereto the slice is copied
            SubConstDFieldRTheta const sub_A = A[ib];
            SubDFieldRTheta sub_A_bar = A_bar[ib];
            DFieldMemRTheta sub_A_alloc(rtheta_mesh, ddc::DeviceAllocator<double>());
            ddc::parallel_deepcopy(sub_A_alloc, sub_A);
            spline_builder(coef, get_const_field(sub_A_alloc));

            // Loop over r, theta
            ddc::parallel_for_each(
                    rtheta_mesh,
                    KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                        IdxR const ir(irtheta);
                        IdxTheta const itheta(irtheta);
                        double const r = ddc::coordinate(ir);
                        double const theta = ddc::coordinate(itheta);

                        // Average over gyro points
                        double tmp = 0.0;
                        for (std::size_t igyro = 0; igyro < nb_gyro_points; igyro++) {
                            // Compute the particle position in (R, Z) coordinate
                            double const alpha = M_PI * 2.0 / static_cast<double>(nb_gyro_points)
                                                 * static_cast<double>(igyro);
                            double const R_p = r * Kokkos::cos(theta)
                                               + rho_L(ir, itheta) * Kokkos::cos(alpha);
                            double const Z_p = r * Kokkos::sin(theta)
                                               + rho_L(ir, itheta) * Kokkos::sin(alpha);

                            // Convert from (R, Z) into (r, theta) coordinate
                            CoordRTheta p = coordinate_transform(R_p, Z_p);

                            // Spline interpolation in (r, theta) coordinate
                            tmp += spline_evaluator(p, get_const_field(coef));
                        }
                        sub_A_bar(ir, itheta) = tmp / static_cast<double>(nb_gyro_points);
                    });

            // Apply periodic boundary condition in theta direction
            ddc::parallel_deepcopy(sub_A_bar[theta_domain.back()], sub_A_bar[theta_domain.front()]);
        });
    }
};

} // namespace detail
