#include <cmath>
#include <random>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/polar_bsplines.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "test_utils.hpp"

template <class T>
struct PolarBsplineFixture;

template <std::size_t D, int C, bool Uniform>
struct PolarBsplineFixture<std::tuple<
        std::integral_constant<std::size_t, D>,
        std::integral_constant<int, C>,
        std::integral_constant<bool, Uniform>>> : public testing::Test
{
    struct DimR
    {
        static constexpr bool PERIODIC = false;
    };
    struct DimP
    {
        static constexpr bool PERIODIC = true;
    };
    struct DimX
    {
        static constexpr bool PERIODIC = false;
    };
    struct DimY
    {
        static constexpr bool PERIODIC = false;
    };
    static constexpr std::size_t spline_degree = D;
    static constexpr int continuity = C;
    using BSplineR = ddc::NonUniformBSplines<DimR, D>;
    using BSplineP = std::
            conditional_t<Uniform, ddc::UniformBSplines<DimP, D>, ddc::NonUniformBSplines<DimP, D>>;

    using GrevillePointsR = ddc::GrevilleInterpolationPoints<
            BSplineR,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    using GrevillePointsP = ddc::GrevilleInterpolationPoints<
            BSplineP,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>;

    using IDimR = typename GrevillePointsR::interpolation_mesh_type;
    using IDimP = typename GrevillePointsP::interpolation_mesh_type;
};

using degrees = std::integer_sequence<std::size_t, 1, 2, 3>;
using continuity = std::integer_sequence<int, -1, 0, 1>;
using is_uniform_types = std::tuple<std::true_type, std::false_type>;

using Cases = tuple_to_types_t<cartesian_product_t<degrees, continuity, is_uniform_types>>;

TYPED_TEST_SUITE(PolarBsplineFixture, Cases);

TYPED_TEST(PolarBsplineFixture, PartitionOfUnity)
{
    int constexpr continuity = TestFixture::continuity;
    using DimR = typename TestFixture::DimR;
    using IDimR = typename TestFixture::IDimR;
    using DVectR = ddc::DiscreteVector<IDimR>;
    using DimP = typename TestFixture::DimP;
    using IDimP = typename TestFixture::IDimP;
    using DVectP = ddc::DiscreteVector<IDimP>;
    using DimX = typename TestFixture::DimX;
    using DimY = typename TestFixture::DimY;
    using PolarCoord = ddc::Coordinate<DimR, DimP>;
    using BSplinesR = typename TestFixture::BSplineR;
    using BSplinesP = typename TestFixture::BSplineP;
    using CircToCart = CircularToCartesian<DimX, DimY, DimR, DimP>;
    using BuilderR = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            IDimR,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::SplineSolver::GINKGO,
            IDimR>;
    using BuilderP = ddc::SplineBuilder<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesP,
            IDimP,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::GINKGO,
            IDimP>;
    using SplineRPBuilder = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            IDimR,
            IDimP,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::GINKGO,
            IDimR,
            IDimP>;
    using SplineRPEvaluator = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            IDimR,
            IDimP,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<DimP>,
            ddc::PeriodicExtrapolationRule<DimP>,
            IDimR,
            IDimP>;
    using DiscreteMapping = DiscreteToCartesian<DimX, DimY, SplineRPBuilder, SplineRPEvaluator>;
    using BSplines = PolarBSplines<BSplinesR, BSplinesP, continuity>;
    using CoordR = ddc::Coordinate<DimR>;
    using CoordP = ddc::Coordinate<DimP>;
    using GrevillePointsR = typename TestFixture::GrevillePointsR;
    using GrevillePointsP = typename TestFixture::GrevillePointsP;

    CoordR constexpr r0(0.);
    CoordR constexpr rN(1.);
    CoordP constexpr p0(0.);
    CoordP constexpr pN(2. * M_PI);
    std::size_t constexpr ncells = 20;

    // 1. Create BSplines
    {
        DVectR constexpr npoints(ncells + 1);
        std::vector<CoordR> breaks(npoints);
        const double dr = (rN - r0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordR(r0 + i * dr);
        }
        ddc::init_discrete_space<BSplinesR>(breaks);
    }
    if constexpr (BSplinesP::is_uniform()) {
        ddc::init_discrete_space<BSplinesP>(p0, pN, ncells);
    } else {
        DVectP constexpr npoints(ncells + 1);
        std::vector<CoordP> breaks(npoints);
        const double dp = (pN - p0) / ncells;
        for (int i(0); i < npoints; ++i) {
            breaks[i] = CoordP(p0 + i * dp);
        }
        ddc::init_discrete_space<BSplinesP>(breaks);
    }

    ddc::init_discrete_space<IDimR>(GrevillePointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(GrevillePointsP::get_sampling());
    ddc::DiscreteDomain<IDimR> interpolation_domain_R(GrevillePointsR::get_domain());
    ddc::DiscreteDomain<IDimP> interpolation_domain_P(GrevillePointsP::get_domain());
    ddc::DiscreteDomain<IDimR, IDimP>
            interpolation_domain(interpolation_domain_R, interpolation_domain_P);

    BuilderR builder_r(interpolation_domain_R);
    BuilderP builder_p(interpolation_domain_P);
    SplineRPBuilder builder_rp(interpolation_domain);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<DimP> p_extrapolation_rule;
    SplineRPEvaluator evaluator_rp(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);

    const CircToCart coord_changer;
    DiscreteMapping const mapping
            = DiscreteMapping::analytical_to_discrete(coord_changer, builder_rp, evaluator_rp);
    ddc::init_discrete_space<BSplines>(mapping, builder_r, builder_p);

    int const n_eval = (BSplinesR::degree() + 1) * (BSplinesP::degree() + 1);
    std::size_t const n_test_points = 100;
    double const dr = (rN - r0) / n_test_points;
    double const dp = (pN - p0) / n_test_points;

    for (std::size_t i(0); i < n_test_points; ++i) {
        for (std::size_t j(0); j < n_test_points; ++j) {
            std::array<double, BSplines::n_singular_basis()> singular_data;
            DSpan1D singular_vals(singular_data.data(), BSplines::n_singular_basis());
            std::array<double, n_eval> data;
            DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesP::degree() + 1);

            PolarCoord const test_point(r0 + i * dr, p0 + j * dp);
            ddc::discrete_space<BSplines>().eval_basis(singular_vals, vals, test_point);
            double total(0.0);
            for (std::size_t k(0); k < BSplines::n_singular_basis(); ++k) {
                total += singular_vals(k);
            }
            for (std::size_t k(0); k < BSplinesR::degree() + 1; ++k) {
                for (std::size_t l(0); l < BSplinesP::degree() + 1; ++l) {
                    total += vals(k, l);
                }
            }
            EXPECT_LE(fabs(total - 1.0), 1.0e-15);
        }
    }
}
