// SPDX-License-Identifier: MIT
#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "gyroaverage_operator.hpp"

namespace {

struct R
{
    static constexpr bool PERIODIC = false;
};

struct Theta
{
    static constexpr bool PERIODIC = true;
};

struct Batch;

// --- Spline definitions
int constexpr BSDegreeR = 3;
int constexpr BSDegreeTheta = 3;

struct BSplinesR : ddc::UniformBSplines<R, BSDegreeR>
{
};
struct BSplinesTheta : ddc::UniformBSplines<Theta, BSDegreeTheta>
{
};

ddc::BoundCond constexpr SplineRBoundary = ddc::BoundCond::GREVILLE;
ddc::BoundCond constexpr SplineThetaBoundary = ddc::BoundCond::PERIODIC;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;

struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};
struct GridBatch : UniformGridBase<Batch>
{
};

using CoordR = ddc::Coordinate<R>;
using CoordTheta = ddc::Coordinate<Theta>;
using CoordRTheta = ddc::Coordinate<R, Theta>;

using IdxR = Idx<GridR>;
using IdxTheta = Idx<GridTheta>;
using IdxBatch = Idx<GridBatch>;
using IdxRTheta = Idx<GridR, GridTheta>;
using IdxRBatch = Idx<GridR, GridBatch>;
using IdxRThetaBatch = Idx<GridR, GridTheta, GridBatch>;

using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeBatch = IdxRange<GridBatch>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
using IdxRangeRBatch = IdxRange<GridR, GridBatch>;
using IdxRangeRThetaBatch = IdxRange<GridR, GridTheta, GridBatch>;

using DFieldMemRTheta = DFieldMem<IdxRangeRTheta>;
using DFieldMemRThetaBatch = DFieldMem<IdxRangeRThetaBatch>;
using DFieldRTheta = DField<IdxRangeRTheta>;
using DFieldRThetaBatch = DField<IdxRangeRThetaBatch>;
using DConstFieldRThetaBatch = ConstField<double, IdxRangeRThetaBatch>;

template <typename Field2DType, typename Field3DType>
void initialise(
        Field2DType const& Rcoord,
        Field2DType const& Zcoord,
        Field2DType const& rho_L,
        Field3DType const& A,
        double const B0,
        double const R0,
        double const gyroradius,
        double const kx,
        double const ky)
{
    IdxRangeRTheta const rtheta_mesh = get_idx_range<GridR, GridTheta>(A);

    // Initialise R, Z and rho_L
    ddc::parallel_for_each(
            rtheta_mesh,
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
                rho_L(ir, itheta) = gyroradius / Kokkos::sqrt(B / B0);
            });

    // Initialise A with Fourier mode
    IdxRangeRThetaBatch const rthetabatch_mesh = get_idx_range<GridR, GridTheta, GridBatch>(A);
    ddc::parallel_for_each(
            rthetabatch_mesh,
            KOKKOS_LAMBDA(IdxRThetaBatch const irthetabatch) {
                IdxR const ir(irthetabatch);
                IdxTheta const itheta(irthetabatch);
                IdxBatch const ibatch(irthetabatch);

                // In non-circular case, it should be computed by
                // (Rcoord(0,Ntheta/2) + Rcoord(0,0)) / 2
                double const R0eff = R0;
                A(ir, itheta, ibatch)
                        = Kokkos::cos(kx * (Rcoord(ir, itheta) - R0eff) + ky * Zcoord(ir, itheta));
            });
}

} // namespace

class GyroAverageCircularParamTests : public ::testing::TestWithParam<double>
{
protected:
    void SetUp() override
    {
        // Get reference gyroradius
        m_gyroradius = GetParam();

        // Start of the domain of interest in the R dimension
        double const r_start = 0.;
        double const r_end = m_r0;

        // Number of discretization points in the R dimension
        std::size_t const nb_r_points = 200;

        // Initialisation of the global domain in R
        ddc::init_discrete_space<BSplinesR>(CoordR(r_start), CoordR(r_end), nb_r_points);
        ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());

        IdxRangeR const r_domain = SplineInterpPointsR::get_domain<GridR>();

        // Start of the domain of interest in the Theta dimension
        double const theta_start = 0.;
        double const theta_end = M_PI * 2.0;

        // Number of discretization points in the Theta dimension
        std::size_t const nb_theta_points = 200;

        // Initialisation of the global domain in Theta
        ddc::init_discrete_space<
                BSplinesTheta>(CoordTheta(theta_start), CoordTheta(theta_end), nb_theta_points);
        ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

        IdxRangeTheta const theta_domain = SplineInterpPointsTheta::get_domain<GridTheta>();

        // Start of the domain of interest in the batch dimension (likely a mixture of
        // phi and vpar)
        double const batch_start = 0.;
        double const batch_end = M_PI * 2.0;

        // Number of discretization points in the batch dimension
        std::size_t const nb_batch_points = 2;

        // Initialisation of the global domain in batch
        IdxRangeBatch const batch_domain
                = ddc::init_discrete_space<GridBatch>(GridBatch::init<GridBatch>(
                        ddc::Coordinate<Batch>(batch_start),
                        ddc::Coordinate<Batch>(batch_end),
                        IdxStep<GridBatch>(nb_batch_points)));

        // Allocate members with move assign operator
        IdxRangeRTheta const rtheta_mesh = IdxRangeRTheta(r_domain, theta_domain);
        IdxRangeRThetaBatch const rthetabatch_mesh
                = IdxRangeRThetaBatch(r_domain, theta_domain, batch_domain);
        m_Rcoord_alloc = DFieldMemRTheta(rtheta_mesh, ddc::DeviceAllocator<double>());
        m_Zcoord_alloc = DFieldMemRTheta(rtheta_mesh, ddc::DeviceAllocator<double>());
        m_rho_L_alloc = DFieldMemRTheta(rtheta_mesh, ddc::DeviceAllocator<double>());

        m_A_alloc = DFieldMemRThetaBatch(rthetabatch_mesh, ddc::DeviceAllocator<double>());
        m_A_bar_alloc = DFieldMemRThetaBatch(rthetabatch_mesh, ddc::DeviceAllocator<double>());

        // Initialise R, Z and rho_L
        DFieldRTheta Rcoord = get_field(m_Rcoord_alloc);
        DFieldRTheta Zcoord = get_field(m_Zcoord_alloc);
        DFieldRTheta rho_L = get_field(m_rho_L_alloc);

        m_kx = m_kperp * 1.0 / Kokkos::sqrt(2.0);
        m_ky = m_kx;
        DFieldRThetaBatch A = get_field(m_A_alloc);
        initialise(Rcoord, Zcoord, rho_L, A, m_B0, m_R0, m_gyroradius, m_kx, m_ky);
    }

protected:
    double m_gyroradius;
    std::size_t const m_nb_gyro_points = 32;
    double const m_B0 = 1.0;
    double const m_R0 = 300.0;
    double const m_r0 = 150.0;
    double const m_kperp = M_PI / 6.0;
    double m_kx;
    double m_ky;

    DFieldMemRTheta m_Rcoord_alloc;
    DFieldMemRTheta m_Zcoord_alloc;
    DFieldMemRTheta m_rho_L_alloc;

    DFieldMemRThetaBatch m_A_alloc;
    DFieldMemRThetaBatch m_A_bar_alloc;
};

struct CartesianToPolar
{
    KOKKOS_INLINE_FUNCTION CoordRTheta operator()(double const x, double const y) const
    {
        double r = Kokkos::hypot(x, y);
        double theta = Kokkos::atan2(y, x);
        if (theta < 0)
            theta += M_PI * 2.0;
        return CoordRTheta(r, theta);
    }
};

TEST_P(GyroAverageCircularParamTests, TestPeriodicity)
{
    DConstFieldRThetaBatch A = get_const_field(this->m_A_alloc);
    DFieldRThetaBatch A_bar = get_field(this->m_A_bar_alloc);

    using GyroAverageOperatorType = detail::GyroAverageOperator<
            Kokkos::DefaultExecutionSpace,
            GridR,
            GridTheta,
            IdxRangeRThetaBatch,
            BSplinesR,
            BSplinesTheta,
            CartesianToPolar>;
    GyroAverageOperatorType
            gyroaverage(get_field(this->m_rho_L_alloc), CartesianToPolar(), this->m_nb_gyro_points);
    gyroaverage(A, A_bar);

    auto h_A_bar_alloc = ddc::create_mirror_and_copy(Kokkos::HostSpace {}, A_bar);
    host_t<DFieldRThetaBatch> h_A_bar = get_field(h_A_bar_alloc);

    IdxRangeTheta const theta_domain = get_idx_range<GridTheta>(A_bar);
    IdxRangeRBatch const rbatch_mesh = get_idx_range<GridR, GridBatch>(A_bar);

    ddc::for_each(rbatch_mesh, [&](IdxRBatch const irbatch) {
        IdxR const ir(irbatch);
        IdxBatch const ibatch(irbatch);
        IdxTheta const theta_front = theta_domain.front();
        IdxTheta const theta_back = theta_domain.back();
        double A_theta_0 = h_A_bar(ir, theta_front, ibatch);
        double A_theta_2PI = h_A_bar(ir, theta_back, ibatch);
        EXPECT_LE(std::abs(A_theta_2PI - A_theta_0), 1e-12);
    });
}

TEST_P(GyroAverageCircularParamTests, TestAnalytical)
{
    DConstFieldRThetaBatch A = get_const_field(this->m_A_alloc);
    DFieldRThetaBatch A_bar = get_field(this->m_A_bar_alloc);

    using GyroAverageOperatorType = detail::GyroAverageOperator<
            Kokkos::DefaultExecutionSpace,
            GridR,
            GridTheta,
            IdxRangeRThetaBatch,
            BSplinesR,
            BSplinesTheta,
            CartesianToPolar>;
    GyroAverageOperatorType
            gyroaverage(get_field(this->m_rho_L_alloc), CartesianToPolar(), this->m_nb_gyro_points);
    gyroaverage(A, A_bar);

    auto h_Rcoord_alloc
            = ddc::create_mirror_and_copy(Kokkos::HostSpace {}, get_field(this->m_Rcoord_alloc));
    host_t<DFieldRTheta> h_Rcoord = get_field(h_Rcoord_alloc);
    auto h_Zcoord_alloc
            = ddc::create_mirror_and_copy(Kokkos::HostSpace {}, get_field(this->m_Zcoord_alloc));
    host_t<DFieldRTheta> h_Zcoord = get_field(h_Zcoord_alloc);
    auto h_rho_L_alloc
            = ddc::create_mirror_and_copy(Kokkos::HostSpace {}, get_field(this->m_rho_L_alloc));
    host_t<DFieldRTheta> h_rho_L = get_field(h_rho_L_alloc);
    auto h_A_bar_alloc = ddc::create_mirror_and_copy(Kokkos::HostSpace {}, A_bar);
    host_t<DFieldRThetaBatch> h_A_bar = get_field(h_A_bar_alloc);

    IdxRangeTheta const theta_domain = get_idx_range<GridTheta>(A_bar);
    IdxRangeRThetaBatch const rthetabatch_mesh = get_idx_range<GridR, GridTheta, GridBatch>(A_bar);

    double const r0 = this->m_r0;
    double const R0 = this->m_R0;
    double const kperp = this->m_kperp;
    double const kx = this->m_kx;
    double const ky = this->m_ky;
    int const nb_gyro_points = this->m_nb_gyro_points;
    ddc::for_each(rthetabatch_mesh, [&](IdxRThetaBatch const irthetabatch) {
        IdxR const ir(irthetabatch);
        IdxTheta const itheta(irthetabatch);
        IdxBatch const ibatch(irthetabatch);
        double const r = ddc::coordinate(ir);
        double const R0eff = R0;
        double const gyroradius = h_rho_L(ir, itheta);
        double const kperprho = kperp * gyroradius;
        double const phase = kx * (h_Rcoord(ir, itheta) - R0eff) + ky * h_Zcoord(ir, itheta);
        double const phik
                = std::atan2(h_Zcoord(ir, itheta), h_Rcoord(ir, itheta) - R0eff) - M_PI / 4;
        std::size_t nb_gyro_points_half = m_nb_gyro_points / 2;
        double const sgn = std::pow(-1.0, static_cast<double>(nb_gyro_points_half));
        double const sgn_tmp = nb_gyro_points % 2 == 0 ? sgn : -sgn;
        double const res_th = std::cyl_bessel_j(0, kperprho) * std::cos(phase);
        double const err_J0_th
                = 2.0
                  * (std::cyl_bessel_j(nb_gyro_points, kperprho) * std::cos(nb_gyro_points * phik)
                             * std::cos(phase) * sgn_tmp
                     + std::cyl_bessel_j(2 * nb_gyro_points, kperprho)
                               * std::cos(2 * nb_gyro_points * phik) * std::cos(phase));
        double const err_J0_num = h_A_bar(ir, itheta, ibatch) - res_th;

        // In the outer region, the particle position can be out of small radius (r), where the values are considered to be zero.
        // Also the gyroaveraged value has the periodicity along theta direction (tested separately above)
        // These conditions are not taken into account in the analytical formula, so we do not test these cases.
        if (r < r0 * 0.8 && itheta < theta_domain.back()) {
            EXPECT_LE(std::abs(err_J0_num - err_J0_th), 5e-2);
        }
    });
}

// Parameterisation over Larmor radius
INSTANTIATE_TEST_SUITE_P(
        GyroAverageCircular,
        GyroAverageCircularParamTests,
        ::testing::Values(
                0.55,
                1.1,
                1.65,
                2.2,
                2.75,
                3.3,
                3.85,
                4.4,
                4.95,
                5.5,
                6.05,
                6.6,
                7.15,
                7.7,
                8.25,
                8.8));
