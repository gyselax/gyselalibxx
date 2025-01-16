#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include <gtest/gtest.h>

#include "gauss_legendre_integration.hpp"
#include "quadrature.hpp"
#include "test_utils.hpp"

namespace {

struct X
{
};

template <class T>
struct GaussLegendreFixture;

static std::vector<std::pair<double, double>> test_domains
        = {{0.0, 1.0}, {1.0, 2.0}, {-0.2, 1.5}, {-1.5, -1.0}};

template <std::size_t ORDER, std::size_t DOMAIN_IDX>
struct GaussLegendreFixture<std::tuple<
        std::integral_constant<std::size_t, ORDER>,
        std::integral_constant<std::size_t, DOMAIN_IDX>>> : public testing::Test
{
    struct GLGrid : NonUniformGridBase<X>
    {
    };
    using GaussLegendreBuilder = GaussLegendre<GLGrid, ORDER>;
    static constexpr std::size_t domain_idx = DOMAIN_IDX;
};

using orders = std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10>;
using domain_idx = std::integer_sequence<std::size_t, 0, 1, 2, 3>;

using Cases = tuple_to_types_t<cartesian_product_t<orders, domain_idx>>;

TYPED_TEST_SUITE(GaussLegendreFixture, Cases);

class fn
{
public:
    constexpr fn(std::size_t n) : m_n(n) {}

    template <class Grid1D>
    KOKKOS_FUNCTION double operator()(Idx<Grid1D> idx) const noexcept
    {
        double x = ddc::coordinate(idx);
        double r = 1.0;
        for (std::size_t i = 0; i < m_n; ++i) {
            r *= x;
        }
        return r;
    }

private:
    std::size_t m_n;
};

static std::map<std::size_t, std::string> type_names;

/// This test integrates polynomials of the form x^p for p <= 2*i-1
/// where i is the order of the GaussLegendre integration method.
///
/// For such polynomials, this quadrature rule is exact (truncation
/// error is exactly zero).
TYPED_TEST(GaussLegendreFixture, Integrate)
{
    using GaussLegendreBuilder = typename TestFixture::GaussLegendreBuilder;
    using GLGrid = typename TestFixture::GLGrid;

    Coord<X> x0(test_domains[TestFixture::domain_idx].first);
    Coord<X> x1(test_domains[TestFixture::domain_idx].second);

    std::stringstream oss;
    oss << std::scientific << std::hexfloat;

    std::cout << "integration at order " << GaussLegendreBuilder::order();
    std::cout << std::endl;

    GaussLegendreBuilder const gl({x0, x1});
    DFieldMem<IdxRange<GLGrid>> gauss_legendre_coefficients
            = gl.template gauss_legendre_coefficients<Kokkos::DefaultExecutionSpace>();
    Quadrature integrate(get_const_field(gauss_legendre_coefficients));

    for (std::size_t p = 0; p <= GaussLegendreBuilder::order(); ++p) {
        fn const f(p);
        double const sol_exact = 1.0 / (p + 1) * (std::pow(x1, p + 1) - std::pow(x0, p + 1));
        double const sol_num = integrate(Kokkos::DefaultExecutionSpace(), f);
        double const err = std::fabs((sol_num - sol_exact) / sol_exact);

        bool ok = true;
        if (sol_num != sol_exact) {
            ok = std::log10(err) < -std::numeric_limits<double>::digits10;
        }

        EXPECT_NEAR(sol_num, sol_exact, 1e-10);

        std::cout << sol_exact << " " << sol_num << std::endl;

        oss.str("");
        oss << " of x^" << std::setw(2) << std::left << p;
        oss << ' ';
        oss << std::fixed << std::setprecision(1) << std::right;
        oss << " on the idx_range [" << std::setw(4) << x0 << ", " << std::setw(4) << x1 << "]";
        oss << std::scientific << std::hexfloat;
        oss << ' ';
        oss << std::setw(25) << std::left << sol_num;
        oss << ' ';
        oss << std::setw(25) << std::left << sol_exact;
        std::string str = oss.str();
        oss.str("");
        oss << std::setw(60) << std::left << str;
        oss << (ok ? "PASSED" : "FAILED");
        std::cout << oss.str() << std::endl;
        std::cout << std::endl;
    }
}

} // namespace
