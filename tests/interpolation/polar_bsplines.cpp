#include <cmath>
#include <random>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "discrete_poloidal_cs_spline_mapping.hpp"
#include "discrete_poloidal_cs_spline_mapping_builder.hpp"
#include "polar_bsplines.hpp"
#include "spline_interpolation.hpp"
#include "test_utils.hpp"
#include "view.hpp"

template <class T>
struct PolarBsplineFixture;

template <std::size_t D, int C, bool Uniform>
struct PolarBsplineFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        std::integral_constant<int, C>,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct R_cov;
    struct Theta_cov;
    struct R
    {
        static bool constexpr PERIODIC = false;
        static bool constexpr IS_COVARIANT = false;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = R_cov;
    };

    struct Theta
    {
        static bool constexpr PERIODIC = true;
        static bool constexpr IS_COVARIANT = false;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = Theta_cov;
    };
    struct R_cov
    {
        static bool constexpr PERIODIC = false;
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = false;
        using Dual = R;
    };

    struct Theta_cov
    {
        static bool constexpr PERIODIC = true;
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = false;
        using Dual = Theta;
    };
    struct X
    {
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = X;
    };
    struct Y
    {
        static bool constexpr IS_COVARIANT = true;
        static bool constexpr IS_CONTRAVARIANT = true;
        using Dual = Y;
    };
    static constexpr std::size_t spline_degree = D;
    static constexpr int continuity = C;
    struct BSplinesR : ddc::NonUniformBSplines<R, D, R::PERIODIC>
    {
    };
    struct BSplinesTheta
        : std::conditional_t<
                  Uniform,
                  ddc::UniformBSplines<Theta, D, Theta::PERIODIC>,
                  ddc::NonUniformBSplines<Theta, D, Theta::PERIODIC>>
    {
    };

    using GrevillePointsR = ddc::GrevilleInterpolationPoints<
            BSplinesR,
            ddc::SplineBuilderClosure::GREVILLE,
            ddc::SplineBuilderClosure::GREVILLE>;
    using GrevillePointsTheta = ddc::GrevilleInterpolationPoints<
            BSplinesTheta,
            ddc::SplineBuilderClosure::PERIODIC,
            ddc::SplineBuilderClosure::PERIODIC>;

    struct GridR : GrevillePointsR::interpolation_discrete_dimension_type
    {
    };
    struct GridTheta : GrevillePointsTheta::interpolation_discrete_dimension_type
    {
    };
    struct BSplines : PolarBSplines<BSplinesR, BSplinesTheta, continuity>
    {
    };

    static constexpr Coord<R> r0 = Coord<R>(0.);
    static constexpr Coord<R> rN = Coord<R>(1.);
    static constexpr Coord<Theta> p0 = Coord<Theta>(0.);
    static constexpr Coord<Theta> pN = Coord<Theta>(2. * M_PI);
    static constexpr std::size_t ncells = 20;

    static void SetUpTestSuite()
    {
        using CircToCart = CircularToCartesian<R, Theta, X, Y>;
        using SplineInterpolatorRTheta_host = SplineInterpolator<
                Kokkos::DefaultHostExecutionSpace,
                IdxRange<BSplinesR, BSplinesTheta>,
                IdxRange<GridR, GridTheta>,
                ExtrapolationRule::Null_Null, // radial extrapolation
                ExtrapolationRule::Periodic, // poloidal extrapolation
                SplineBoundaryClosure::Greville_Greville, // radial closure condition
                SplineBoundaryClosure::Periodic>;

        // 1. Create BSplines
        {
            IdxStep<GridR> constexpr npoints(ncells + 1);
            std::vector<Coord<R>> breaks(npoints);
            const double dr = (rN - r0) / ncells;
            for (int i(0); i < npoints; ++i) {
                breaks[i] = Coord<R>(r0 + i * dr);
            }
            ddc::init_discrete_space<BSplinesR>(breaks);
        }
        if constexpr (BSplinesTheta::is_uniform()) {
            ddc::init_discrete_space<BSplinesTheta>(p0, pN, ncells);
        } else {
            IdxStep<GridTheta> constexpr npoints(ncells + 1);
            std::vector<Coord<Theta>> breaks(npoints);
            const double dp = (pN - p0) / ncells;
            for (int i(0); i < npoints; ++i) {
                breaks[i] = Coord<Theta>(p0 + i * dp);
            }
            ddc::init_discrete_space<BSplinesTheta>(breaks);
        }

        ddc::init_discrete_space<GridR>(GrevillePointsR::template get_sampling<GridR>());
        ddc::init_discrete_space<GridTheta>(
                GrevillePointsTheta::template get_sampling<GridTheta>());
        IdxRange<GridR> interpolation_idx_range_r(GrevillePointsR::template get_domain<GridR>());
        IdxRange<GridTheta> interpolation_idx_range_theta(
                GrevillePointsTheta::template get_domain<GridTheta>());
        IdxRange<GridR, GridTheta>
                interpolation_idx_range(interpolation_idx_range_r, interpolation_idx_range_theta);

        SplineInterpolatorRTheta_host interpolator(interpolation_idx_range);

        const CircToCart coord_changer;
        DiscretePoloidalCSSplineMappingBuilder<X, Y, SplineInterpolatorRTheta_host>
                mapping_builder(Kokkos::DefaultHostExecutionSpace(), coord_changer, interpolator);
        DiscretePoloidalCSSplineMapping mapping = mapping_builder();
        ddc::init_discrete_space<BSplines>(mapping);
    }
};

using degrees = std::integer_sequence<std::size_t, 1, 2, 3>;
using continuity = std::integer_sequence<int, -1, 0, 1>;
using is_uniform_types = std::tuple<std::true_type, std::false_type>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees, continuity, is_uniform_types>>;

TYPED_TEST_SUITE(PolarBsplineFixture, Cases);

TYPED_TEST(PolarBsplineFixture, PartitionOfUnity)
{
    using R = typename TestFixture::R;
    using Theta = typename TestFixture::Theta;
    using PolarCoord = Coord<R, Theta>;
    using BSplinesR = typename TestFixture::BSplinesR;
    using BSplinesTheta = typename TestFixture::BSplinesTheta;
    using BSplines = typename TestFixture::BSplines;

    Coord<R> const r0 = TestFixture::r0;
    Coord<R> const rN = TestFixture::rN;
    Coord<Theta> const p0 = TestFixture::p0;
    Coord<Theta> const pN = TestFixture::pN;

    int const n_eval = (BSplinesR::degree() + 1) * (BSplinesTheta::degree() + 1);
    std::size_t const n_test_points = 100;
    double const dr = (rN - r0) / n_test_points;
    double const dp = (pN - p0) / n_test_points;

    for (std::size_t i(0); i < n_test_points; ++i) {
        for (std::size_t j(0); j < n_test_points; ++j) {
            std::array<double, BSplines::n_singular_basis()> singular_data;
            DSpan1D singular_vals(singular_data.data(), BSplines::n_singular_basis());
            std::array<double, n_eval> data;
            DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesTheta::degree() + 1);

            PolarCoord const test_point(r0 + i * dr, p0 + j * dp);
            ddc::discrete_space<BSplines>().eval_basis(singular_vals, vals, test_point);
            double total(0.0);
            for (std::size_t k(0); k < BSplines::n_singular_basis(); ++k) {
                total += singular_vals(k);
            }
            for (std::size_t k(0); k < BSplinesR::degree() + 1; ++k) {
                for (std::size_t l(0); l < BSplinesTheta::degree() + 1; ++l) {
                    total += vals(k, l);
                }
            }
            EXPECT_LE(fabs(total - 1.0), 1.0e-15);
        }
    }
}

