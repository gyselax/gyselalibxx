#include <cmath>
#include <random>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_builder.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>
#include <sll/test_utils.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "polar_bsplines.hpp"

template <class T>
struct PolarBsplineFixture;

template <std::size_t D, int C, bool Uniform>
struct PolarBsplineFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        std::integral_constant<int, C>,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct R
    {
        static constexpr bool PERIODIC = false;
    };
    struct Theta
    {
        static constexpr bool PERIODIC = true;
    };
    struct X
    {
        static constexpr bool PERIODIC = false;
    };
    struct Y
    {
        static constexpr bool PERIODIC = false;
    };
    static constexpr std::size_t spline_degree = D;
    static constexpr int continuity = C;
    struct BSplinesR : ddc::NonUniformBSplines<R, D>
    {
    };
    struct BSplinesTheta
        : std::conditional_t<
                  Uniform,
                  ddc::UniformBSplines<Theta, D>,
                  ddc::NonUniformBSplines<Theta, D>>
    {
    };

    using GrevillePointsR = ddc::GrevilleInterpolationPoints<
            BSplinesR,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    using GrevillePointsTheta = ddc::GrevilleInterpolationPoints<
            BSplinesTheta,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>;

    struct GridR : GrevillePointsR::interpolation_discrete_dimension_type
    {
    };
    struct GridTheta : GrevillePointsTheta::interpolation_discrete_dimension_type
    {
    };
    struct BSplines : PolarBSplines<BSplinesR, BSplinesTheta, continuity>
    {
    };
};

using degrees = std::integer_sequence<std::size_t, 1, 2, 3>;
using continuity = std::integer_sequence<int, -1, 0, 1>;
using is_uniform_types = std::tuple<std::true_type, std::false_type>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees, continuity, is_uniform_types>>;

TYPED_TEST_SUITE(PolarBsplineFixture, Cases);

TYPED_TEST(PolarBsplineFixture, PartitionOfUnity)
{
    using R = typename TestFixture::R;
    using GridR = typename TestFixture::GridR;
    using IdxStepR = IdxStep<GridR>;
    using Theta = typename TestFixture::Theta;
    using GridTheta = typename TestFixture::GridTheta;
    using IdxStepTheta = IdxStep<GridTheta>;
    using X = typename TestFixture::X;
    using Y = typename TestFixture::Y;
    using PolarCoord = Coord<R, Theta>;
    using BSplinesR = typename TestFixture::BSplinesR;
    using BSplinesTheta = typename TestFixture::BSplinesTheta;
    using CircToCart = CircularToCartesian<R, Theta, X, Y>;
    using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridTheta>;
    using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;
    using BSplines = typename TestFixture::BSplines;
    using CoordR = Coord<R>;
    using CoordTheta = Coord<Theta>;
    using GrevillePointsR = typename TestFixture::GrevillePointsR;
    using GrevillePointsTheta = typename TestFixture::GrevillePointsTheta;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordTheta constexpr p0(0.);
    CoordTheta constexpr pN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        IdxStepR constexpr npoints(ncells + 1);
        std::vector<CoordR> breaks(npoints);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks);
    }
    if constexpr (BSplinesTheta::is_uniform()) {
        ddc::init_discrete_space<BSplinesTheta>(p0, pN, ncells);
    } else {
        IdxStepTheta constexpr npoints(ncells + 1);
        std::vector<CoordTheta> breaks(npoints);
        const double dp = (pN - p0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordTheta(p0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesTheta>(breaks);
    }

    ddc::init_discrete_space<GridR>(GrevillePointsR::template get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(GrevillePointsTheta::template get_sampling<GridTheta>());
    IdxRange<GridR> interpolation_idx_range_r(GrevillePointsR::template get_domain<GridR>());
    IdxRange<GridTheta> interpolation_idx_range_theta(
            GrevillePointsTheta::template get_domain<GridTheta>());
    IdxRange<GridR, GridTheta>
            interpolation_idx_range(interpolation_idx_range_r, interpolation_idx_range_theta);

    SplineRThetaBuilder_host builder_rp(interpolation_idx_range);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator evaluator_rp(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    const CircToCart coord_changer;
    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator>
            mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    coord_changer,
                    builder_rp,
                    evaluator_rp);
    DiscreteToCartesian mapping = mapping_builder();
    ddc::init_discrete_space<BSplines>(mapping);

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
