// SPDX-License-Identifier: MIT
#include <cmath>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "gyroaverage_operator.hpp"
#include "spline_interpolation.hpp"
#include "test_math_utils.hpp"



namespace {

struct R_cov;
struct Theta_cov;

struct R
{
    static constexpr bool PERIODIC = false;
    /// The corresponding type in the dual space.
    using Dual = R_cov;
};

struct Theta
{
    static constexpr bool PERIODIC = true;
    /// The corresponding type in the dual space.
    using Dual = Theta_cov;
};

struct Mu;

/// A dummy dimension used to check batching
struct Batch;

// --- Spline definitions
int constexpr BSDegreeR = 3;
int constexpr BSDegreeTheta = 3;

struct BSplinesR : ddc::UniformBSplines<R, BSDegreeR, R::PERIODIC>
{
};
struct BSplinesTheta : ddc::UniformBSplines<Theta, BSDegreeTheta, Theta::PERIODIC>
{
};

ddc::SplineBuilderClosure constexpr SplineRClosure = ddc::SplineBuilderClosure::GREVILLE;
ddc::SplineBuilderClosure constexpr SplineThetaClosure = ddc::SplineBuilderClosure::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRClosure, SplineRClosure>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaClosure, SplineThetaClosure>;

struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};
struct GridMu : NonUniformGridBase<Mu>
{
};
struct GridBatch : UniformGridBase<Batch>
{
};

template <class ExecutionSpace>
using SplineRThetaInterpolatorType = SplineInterpolator<
        ExecutionSpace,
        IdxRange<BSplinesR, BSplinesTheta>,
        IdxRange<GridR, GridTheta>,
        ExtrapolationRule::Null_Null, // radial extrapolation
        ExtrapolationRule::Periodic, // poloidal extrapolation
        SplineBoundaryClosures<
                SplineRClosure, // boundary at r=0
                SplineRClosure>, // boundary at rmax
        SplineBoundaryClosures<SplineThetaClosure, SplineThetaClosure>>;

using CoordR = Coord<R>;
using CoordTheta = Coord<Theta>;
using CoordRTheta = Coord<R, Theta>;
using CoordMu = Coord<Mu>;

using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxMu = Idx<GridMu>;
using IdxRTheta = Idx<GridR, GridTheta>;
using IdxMuR = Idx<GridMu, GridR>;
using IdxMuRTheta = Idx<GridMu, GridR, GridTheta>;
using IdxBatch = Idx<GridBatch>;
using IdxMuBatchRTheta = Idx<GridMu, GridBatch, GridR, GridTheta>;

using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeMu = IdxRange<GridMu>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
using IdxRangeMuR = IdxRange<GridMu, GridR>;
using IdxRangeMuRTheta = IdxRange<GridMu, GridR, GridTheta>;
using IdxRangeBatch = IdxRange<GridBatch>;
using IdxRangeMuBatchRTheta = IdxRange<GridMu, GridBatch, GridR, GridTheta>;

using DFieldMemMu = DFieldMem<IdxRangeMu>;
using DFieldMemRTheta = DFieldMem<IdxRangeRTheta>;
using DFieldMemMuRTheta = DFieldMem<IdxRangeMuRTheta>;
using DFieldMemMuBatchRTheta = DFieldMem<IdxRangeMuBatchRTheta>;
using DFieldMu = DField<IdxRangeMu>;
using DFieldRTheta = DField<IdxRangeRTheta>;
using DFieldMuRTheta = DField<IdxRangeMuRTheta>;
using DFieldMuBatchRTheta = DField<IdxRangeMuBatchRTheta>;
using DConstFieldRTheta = ConstField<double, IdxRangeRTheta>;

void initialise(
        DFieldRTheta const& Rcoord,
        DFieldRTheta const& Zcoord,
        DFieldRTheta const& larmor_factor_space,
        DFieldMu const& larmor_factor_params,
        DFieldRTheta const& A,
        double const B0,
        double const R0,
        double const kx,
        double const ky)
{
    IdxRangeRTheta const rtheta_idx_range = get_idx_range(larmor_factor_space);
    IdxRangeMu const mu_idx_range = get_idx_range(larmor_factor_params);

    // Initialise the mu part of the Larmor radius
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            mu_idx_range,
            KOKKOS_LAMBDA(IdxMu const imu) {
                double const mu = ddc::coordinate(imu);
                // Single species with mass and charge equal to 1
                double const mass = 1.0;
                double const charge = 1.0;
                // rho = sqrt(2 m_s mu / B) / q_s (space dependence through B handled below)
                larmor_factor_params(imu) = Kokkos::sqrt(2 * mu * mass) / Kokkos::fabs(charge);
            });

    // Initialise R, Z and the spatial part of the Larmor radius
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            rtheta_idx_range,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                IdxR const ir(irtheta);
                IdxTheta const itheta(irtheta);
                double const r = ddc::coordinate(ir);
                double const theta = ddc::coordinate(itheta);
                double const R = R0 + r * Kokkos::cos(theta);
                double const Z = r * Kokkos::sin(theta);
                double const B = B0 * R0 / R;
                Rcoord(ir, itheta) = R;
                Zcoord(ir, itheta) = Z;
                larmor_factor_space(ir, itheta) = Kokkos::sqrt(B0 / B);
            });

    // Initialise A with Fourier mode
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            rtheta_idx_range,
            KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                IdxR const ir(irtheta);
                IdxTheta const itheta(irtheta);

                // In non-circular case, R0 effective should be computed as
                // R0eff = (Rcoord(0,Ntheta/2) + Rcoord(0,0)) / 2
                double const R0eff = R0;
                A(ir, itheta)
                        = Kokkos::cos(kx * (Rcoord(ir, itheta) - R0eff) + ky * Zcoord(ir, itheta));
            });
}

} // namespace

/**
 * @brief Evaluate the accuracy of GyroAverageOperator on a polar mesh.
 *
 * The test considers a single species (not indexed) with m_s = q_s = 1 a multiple values of mu.
 * The analytical formula used for validation takes into account the number of sample points and
 * is an extension of the one in [K. Obrejan, PFR 10 (2015)]. The remaining error mainly comes from
 * the interpolation at the sampling points.
 */
class GyroAverageCircularParamTests : public ::testing::TestWithParam<std::vector<double>>
{
protected:
    void SetUp() override
    {
        // Get the values of the magnetic moment
        std::vector<double> const mu_values = GetParam();

        // Start of the domain of interest in the R dimension
        CoordR const r_start(0.);
        CoordR const r_end(m_r0);

        // Number of discretization points in the R dimension
        std::size_t const nb_r_points = 128;

        // Initialisation of the global domain in R
        ddc::init_discrete_space<BSplinesR>(r_start, r_end, nb_r_points);
        ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());

        IdxRangeR const r_domain = SplineInterpPointsR::get_domain<GridR>();

        // Start of the domain of interest in the Theta dimension
        CoordTheta const theta_start(0.);
        CoordTheta const theta_end(M_PI * 2.0);

        // Number of discretization points in the Theta dimension
        std::size_t const nb_theta_points = 512;

        // Initialisation of the global domain in Theta
        ddc::init_discrete_space<BSplinesTheta>(theta_start, theta_end, nb_theta_points);
        ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

        IdxRangeTheta const theta_domain = SplineInterpPointsTheta::get_domain<GridTheta>();

        // Initialisation of the global domain in mu
        std::vector<CoordMu> mu_coords(mu_values.begin(), mu_values.end());
        ddc::init_discrete_space<GridMu>(GridMu::init<GridMu>(mu_coords));
        IdxRangeMu const mu_domain(IdxMu(0), IdxStep<GridMu>(mu_coords.size()));

        // Allocate members with move assign operator
        IdxRangeRTheta const rtheta_mesh = IdxRangeRTheta(r_domain, theta_domain);
        IdxRangeMuRTheta const rthetamu_mesh = IdxRangeMuRTheta(r_domain, theta_domain, mu_domain);
        m_Rcoord_alloc = DFieldMemRTheta(rtheta_mesh);
        m_Zcoord_alloc = DFieldMemRTheta(rtheta_mesh);
        m_larmor_factor_space_alloc = DFieldMemRTheta(rtheta_mesh);
        m_larmor_factor_params_alloc = DFieldMemMu(mu_domain);

        m_A_alloc = DFieldMemRTheta(rtheta_mesh);
        m_A_gavg_alloc = DFieldMemMuRTheta(rthetamu_mesh);

        // Corresponding fields
        DFieldRTheta Rcoord = get_field(m_Rcoord_alloc);
        DFieldRTheta Zcoord = get_field(m_Zcoord_alloc);
        DFieldRTheta larmor_factor_space = get_field(m_larmor_factor_space_alloc);
        DFieldMu larmor_factor_params = get_field(m_larmor_factor_params_alloc);
        DFieldRTheta A = get_field(m_A_alloc);

        // Initialise data using kernels on the device
        m_kx = m_kperp * 1.0 / Kokkos::sqrt(2.0);
        m_ky = m_kx;
        initialise(
                Rcoord,
                Zcoord,
                larmor_factor_space,
                larmor_factor_params,
                A,
                m_B0,
                m_R0,
                m_kx,
                m_ky);
    }

    /**
     * @brief Check the values stored in m_A_gavg_alloc against the analytical gyroaverage of A.
     */
    void check_against_analytical();

protected:
    std::size_t const m_nb_gyro_points = 6; // Must be even for the error estimate below to apply
    double const m_B0 = 1.0;
    double const m_R0 = 300.0;
    double const m_r0 = 150.0;
    double const m_kperp = M_PI / 6.0;
    double m_kx = 0.0;
    double m_ky = 0.0;

    DFieldMemRTheta m_Rcoord_alloc;
    DFieldMemRTheta m_Zcoord_alloc;
    DFieldMemRTheta m_larmor_factor_space_alloc;
    DFieldMemMu m_larmor_factor_params_alloc;

    DFieldMemRTheta m_A_alloc;
    DFieldMemMuRTheta m_A_gavg_alloc;
};

using CartesianToPolar = CartesianToCircular<X_pC, Y_pC, R, Theta>;

using SplineRThetaInterpolator = SplineRThetaInterpolatorType<Kokkos::DefaultExecutionSpace>;

using GyroAverageOperatorType = GyroAverageOperator<SplineRThetaInterpolator, IdxRangeMu>;

using GyroAverageOperatorBatchedType
        = GyroAverageOperator<SplineRThetaInterpolator, IdxRangeMu, IdxRangeBatch>;

TEST_P(GyroAverageCircularParamTests, TestAnalytical)
{
    DFieldRTheta A = get_field(m_A_alloc);
    DFieldMuRTheta A_gavg = get_field(m_A_gavg_alloc);

    IdxRangeMuRTheta const rthetamu_idx_range = get_idx_range(A_gavg);
    IdxRangeRTheta const rtheta_idx_range(rthetamu_idx_range);
    SplineRThetaInterpolator const spline_interpolator(rtheta_idx_range);

    auto [gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc]
            = GyroAverageOperatorType::build_circular_coefficients(
                    CartesianToPolar(),
                    get_field(m_larmor_factor_params_alloc),
                    get_field(m_larmor_factor_space_alloc),
                    m_nb_gyro_points);

    GyroAverageOperatorType
            gyroaverage(spline_interpolator, gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc);

    // Apply the gyroaverage for each value of mu
    IdxRangeMu const mu_idx_range = get_idx_range<GridMu>(A_gavg);
    ddc::host_for_each(mu_idx_range, [&](IdxMu const imu) {
        ddc::parallel_deepcopy(A_gavg[imu], A);
        gyroaverage(A_gavg[imu], imu);
    });

    Kokkos::fence();

    check_against_analytical();
}

TEST_P(GyroAverageCircularParamTests, TestAnalyticalBatched)
{
    IdxRangeBatch const batch_idx_range(IdxBatch(0), IdxStep<GridBatch>(2));

    DFieldRTheta A = get_field(m_A_alloc);
    DFieldMuRTheta A_gavg = get_field(m_A_gavg_alloc);

    IdxRangeMuRTheta const rthetamu_idx_range = get_idx_range(A_gavg);
    IdxRangeRTheta const rtheta_idx_range(rthetamu_idx_range);
    SplineRThetaInterpolator const spline_interpolator(rtheta_idx_range);

    auto [gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc]
            = GyroAverageOperatorBatchedType::build_circular_coefficients(
                    CartesianToPolar(),
                    get_field(m_larmor_factor_params_alloc),
                    get_field(m_larmor_factor_space_alloc),
                    m_nb_gyro_points);

    GyroAverageOperatorBatchedType
            gyroaverage(spline_interpolator, gyrocoeffs_coords_alloc, gyrocoeffs_weights_alloc);

    // Copy A identically to every (mu, batch) slice
    DFieldMemMuBatchRTheta A_gavg_batched_alloc(
            IdxRangeMuBatchRTheta(batch_idx_range, rthetamu_idx_range));
    DFieldMuBatchRTheta A_gavg_batched = get_field(A_gavg_batched_alloc);
    ddc::parallel_copy(Kokkos::DefaultExecutionSpace(), A_gavg_batched, get_const_field(A));

    // Apply the gyroaverage to all values of mu and batch at once
    gyroaverage(A_gavg_batched);

    // Each batch slice should match the analytical result
    ddc::host_for_each(batch_idx_range, [&](IdxBatch const ibatch) {
        ddc::parallel_deepcopy(A_gavg, A_gavg_batched[ibatch]);
        Kokkos::fence();
        check_against_analytical();
    });
}

void GyroAverageCircularParamTests::check_against_analytical()
{
    DFieldMuRTheta A_gavg = get_field(m_A_gavg_alloc);
    IdxRangeMuRTheta const rthetamu_idx_range = get_idx_range(A_gavg);
    IdxRangeMu const mu_idx_range(rthetamu_idx_range);

    auto Rcoord_alloc_host = ddc::create_mirror_and_copy(get_field(m_Rcoord_alloc));
    host_t<DFieldRTheta> Rcoord_host = get_field(Rcoord_alloc_host);
    auto Zcoord_alloc_host = ddc::create_mirror_and_copy(get_field(m_Zcoord_alloc));
    host_t<DFieldRTheta> Zcoord_host = get_field(Zcoord_alloc_host);
    auto larmor_factor_space_alloc_host
            = ddc::create_mirror_and_copy(get_field(m_larmor_factor_space_alloc));
    host_t<DFieldRTheta> larmor_factor_space_host = get_field(larmor_factor_space_alloc_host);
    auto larmor_factor_params_alloc_host
            = ddc::create_mirror_and_copy(get_field(m_larmor_factor_params_alloc));
    host_t<DFieldMu> larmor_factor_params_host = get_field(larmor_factor_params_alloc_host);
    auto A_gavg_alloc_host = ddc::create_mirror_and_copy(A_gavg);
    host_t<DFieldMuRTheta> A_gavg_host = get_field(A_gavg_alloc_host);

    IdxRangeR const r_idx_range(rthetamu_idx_range);
    IdxRangeTheta const theta_idx_range(rthetamu_idx_range);

    IdxStep<GridMu, GridR, GridTheta> rthetamu_check_strides(1, 5, 11);
    IdxStep<GridMu, GridR, GridTheta> rthetamu_check_extents(
            mu_idx_range.size(),
            r_idx_range.size() / 5,
            theta_idx_range.size() / 11);
    IdxRangeSlice<GridMu, GridR, GridTheta> rthetamu_idx_range_check(
            rthetamu_idx_range.front(),
            rthetamu_check_extents,
            rthetamu_check_strides);

    double const r0 = m_r0;
    double const R0 = m_R0;
    double const kperp = m_kperp;
    double const kx = m_kx;
    double const ky = m_ky;
    int const nb_gyro_points = m_nb_gyro_points;
    ddc::host_for_each(rthetamu_idx_range_check, [&](IdxMuRTheta const irthetamu) {
        IdxR const ir(irthetamu);
        IdxTheta const itheta(irthetamu);
        IdxMu const imu(irthetamu);
        double const r = ddc::coordinate(ir);
        // In the outer region, the sample points will be outside of the domain, where the boundary
        // conditions depend on the interpolator.
        // This is not taken into account in the analytical formula, so we do not test these cases.
        if (r < r0 * 0.8) {
            double const R0eff = R0;
            double const gyroradius
                    = larmor_factor_params_host(imu) * larmor_factor_space_host(ir, itheta);
            double const kperprho = kperp * gyroradius;
            double const phase
                    = kx * (Rcoord_host(ir, itheta) - R0eff) + ky * Zcoord_host(ir, itheta);
            double const res_th = cyl_bessel_j(0, kperprho) * std::cos(phase);
            double const err_J0_num = A_gavg_host(ir, itheta, imu) - res_th;

            double const phik = M_PI / 4;
            std::size_t nb_gyro_points_half = m_nb_gyro_points / 2;
            double const sgn = std::pow(-1.0, static_cast<double>(nb_gyro_points_half));
            double const err_J0_th = 2.0 * std::cos(phase)
                                     * (cyl_bessel_j(nb_gyro_points, kperprho)
                                                * std::cos(nb_gyro_points * phik) * sgn
                                        + cyl_bessel_j(2 * nb_gyro_points, kperprho)
                                                  * std::cos(2 * nb_gyro_points * phik));
            EXPECT_NEAR(err_J0_num, err_J0_th, 1e-3);
        }
    });
}

INSTANTIATE_TEST_SUITE_P(
        GyroAverageCircular,
        GyroAverageCircularParamTests,
        ::testing::Values(std::vector<double> {0.4, 10.0, 30.0}));
