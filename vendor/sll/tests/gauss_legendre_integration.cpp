#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>

#include <sll/gauss_legendre_integration.hpp>

#include <gtest/gtest.h>

class fn
{
public:
    constexpr fn(std::size_t n) : m_n(n) {}

    constexpr double operator()(double x) const noexcept
    {
        double r = 1.0;
        for (std::size_t i = 0; i < m_n; ++i) {
            r *= x;
        }
        return r;
    }

private:
    std::size_t m_n;
};

struct X
{
};

static std::map<std::size_t, std::string> type_names;

/// This test integrates polynomials of the form x^p for p <= 2*i-1
/// where i is the order of the GaussLegendre integration method.
///
/// For such polynomials, this quadrature rule is exact (truncation
/// error is exactly zero).
int test_integrate()
{
    std::stringstream oss;
    oss << std::scientific << std::hexfloat;

    std::vector<std::pair<double, double>> idx_ranges
            = {{0.0, 1.0}, {1.0, 2.0}, {-0.2, 1.5}, {-1.5, -1.0}};

    bool test_passed = true;

    for (std::size_t order = 1; order <= GaussLegendre<X>::max_order(); ++order) {
        GaussLegendre<X> const gl(order);

        std::cout << "integration at order " << order;
        std::cout << std::endl;

        for (std::size_t p = 0; p < 2 * order; ++p) {
            fn const f(p);
            for (std::size_t i = 0; i < idx_ranges.size(); ++i) {
                double const sol_exact = 1.0 / (p + 1)
                                         * (std::pow(idx_ranges[i].second, p + 1)
                                            - std::pow(idx_ranges[i].first, p + 1));
                double const sol_num = gl.integrate(f, idx_ranges[i].first, idx_ranges[i].second);
                double const err = std::fabs((sol_num - sol_exact) / sol_exact);

                bool ok = true;
                if (sol_num != sol_exact) {
                    ok = std::log10(err) < -std::numeric_limits<double>::digits10;
                }

                test_passed = test_passed && ok;

                oss.str("");
                oss << " of x^" << std::setw(2) << std::left << p;
                oss << ' ';
                oss << std::fixed << std::setprecision(1) << std::right;
                oss << " on the idx_range [" << std::setw(4) << idx_ranges[i].first << ", "
                    << std::setw(4) << idx_ranges[i].second << "]";
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
            }
        }
        std::cout << std::endl;
    }

    return test_passed;
}

/// This test integrates polynomials of the form x^p for p <= 2*order-1
/// where order is the order of the GaussLegendre integration method.
///
/// For such polynomials, this quadrature rule is exact (truncation
/// error is exactly zero).
int test_compute_points_and_weights()
{
    std::stringstream oss;
    oss << std::scientific << std::hexfloat;

    std::vector<std::pair<double, double>> idx_ranges
            = {{0.0, 1.0}, {1.0, 2.0}, {-0.2, 1.5}, {-1.5, -1.0}};

    bool test_passed = true;

    struct GridX
    {
    };

    for (std::size_t order = 1; order <= GaussLegendre<X>::max_order(); ++order) {
        ddc::DiscreteDomain<GridX>
                idx_range(ddc::DiscreteElement<GridX> {0}, ddc::DiscreteVector<GridX> {order});
        ddc::Chunk<ddc::Coordinate<X>, ddc::DiscreteDomain<GridX>> gl_points(idx_range);
        ddc::Chunk<double, ddc::DiscreteDomain<GridX>> gl_weights(idx_range);
        GaussLegendre<X> const gl(order);

        std::cout << "integration at order " << order;
        std::cout << std::endl;

        for (std::size_t p = 0; p < 2 * order; ++p) {
            fn const f(p);
            for (std::size_t i = 0; i < idx_ranges.size(); ++i) {
                gl.compute_points_and_weights(
                        gl_points.span_view(),
                        gl_weights.span_view(),
                        ddc::Coordinate<X>(idx_ranges[i].first),
                        ddc::Coordinate<X>(idx_ranges[i].second));
                double const sol_exact = 1.0 / (p + 1)
                                         * (std::pow(idx_ranges[i].second, p + 1)
                                            - std::pow(idx_ranges[i].first, p + 1));
                double sol_num = 0.0;
                for (auto xi : idx_range) {
                    sol_num += gl_weights(xi) * f(gl_points(xi));
                }
                double const err = std::fabs((sol_num - sol_exact) / sol_exact);

                bool ok = true;
                if (sol_num != sol_exact) {
                    ok = std::log10(err) < -std::numeric_limits<double>::digits10;
                }

                test_passed = test_passed && ok;

                oss.str("");
                oss << " of x^" << std::setw(2) << std::left << p;
                oss << ' ';
                oss << std::fixed << std::setprecision(1) << std::right;
                oss << " on the idx_range [" << std::setw(4) << idx_ranges[i].first << ", "
                    << std::setw(4) << idx_ranges[i].second << "]";
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
            }
        }
        std::cout << std::endl;
    }

    return test_passed;
}

TEST(GaussLegendre, IntegrateDouble)
{
    ASSERT_TRUE(test_integrate());
    ASSERT_TRUE(test_compute_points_and_weights());
}
