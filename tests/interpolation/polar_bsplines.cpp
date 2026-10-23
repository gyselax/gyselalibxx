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

TYPED_TEST(PolarBsplineFixture, BasisAndNDerivs)
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

    std::size_t constexpr nr = BSplinesR::degree() + 1;
    std::size_t constexpr ntheta = BSplinesTheta::degree() + 1;
    std::size_t constexpr n_singular = BSplines::n_singular_basis();
    int constexpr n_derivs_r = BSplinesR::degree();
    int constexpr n_derivs_theta = BSplinesTheta::degree();

    std::size_t const n_test_points = 20;
    double const dr = (rN - r0) / n_test_points;
    double const dp = (pN - p0) / n_test_points;

    for (std::size_t i(0); i < n_test_points; ++i) {
        for (std::size_t j(0); j < n_test_points; ++j) {
            PolarCoord const test_point(r0 + i * dr, p0 + j * dp);

            std::array<double, n_singular * nr * ntheta> singular_data;
            DSpan3D singular_vals_and_derivs(
                    singular_data.data(),
                    n_singular,
                    n_derivs_r + 1,
                    n_derivs_theta + 1);
            std::array<double, nr * ntheta * nr * ntheta> data;
            DSpan4D vals_and_derivs(data.data(), nr, ntheta, n_derivs_r + 1, n_derivs_theta + 1);
            Idx<BSplinesR, BSplinesTheta> idx_front
                    = ddc::discrete_space<BSplines>().eval_basis_and_n_derivs(
                            singular_vals_and_derivs,
                            vals_and_derivs,
                            test_point,
                            n_derivs_r,
                            n_derivs_theta);

            for (int a(0); a <= n_derivs_r; ++a) {
                for (int b(0); b <= n_derivs_theta; ++b) {
                    std::array<double, n_singular> singular_ref_data;
                    DSpan1D singular_ref(singular_ref_data.data(), n_singular);
                    std::array<double, nr * ntheta> ref_data;
                    DSpan2D ref(ref_data.data(), nr, ntheta);
                    Idx<BSplinesR, BSplinesTheta> idx_front_ref;
                    if (a == 0 && b == 0) {
                        idx_front_ref = ddc::discrete_space<BSplines>()
                                                .eval_basis(singular_ref, ref, test_point);
                    } else {
                        idx_front_ref = ddc::discrete_space<BSplines>().eval_deriv(
                                singular_ref,
                                ref,
                                test_point,
                                Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(a, b));
                    }
                    EXPECT_EQ(idx_front, idx_front_ref);
                    for (std::size_t k(0); k < n_singular; ++k) {
                        EXPECT_NEAR(
                                singular_vals_and_derivs(k, a, b),
                                singular_ref(k),
                                1e-12 * (1.0 + std::fabs(singular_ref(k))));
                    }
                    for (std::size_t k(0); k < nr; ++k) {
                        for (std::size_t l(0); l < ntheta; ++l) {
                            EXPECT_NEAR(
                                    vals_and_derivs(k, l, a, b),
                                    ref(k, l),
                                    1e-12 * (1.0 + std::fabs(ref(k, l))));
                        }
                    }
                }
            }
        }
    }
}
