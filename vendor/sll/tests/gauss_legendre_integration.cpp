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

    bool test_passed = true;

    for (std::size_t i = 1; i <= GaussLegendre::max_order(); ++i) {
        GaussLegendre const gl(i);

        std::cout << "integration at order " << i;
        std::cout << std::endl;

        for (std::size_t p = 0; p < 2 * i; ++p) {
            fn const f(p);
            double const sol_exact = 1.0 / (p + 1);
            double const sol_num = gl.integrate(f, 0, 1);
            double const err = std::fabs(sol_num - sol_exact) / sol_exact;

            bool ok = true;
            if (sol_num != sol_exact) {
                ok = std::log10(err) < -std::numeric_limits<double>::digits10;
            }

            test_passed = test_passed && ok;

            oss.str("");
            oss << " of x^" << std::setw(2) << std::left << p;
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
        std::cout << std::endl;
    }

    return test_passed;
};

TEST(GaussLegendre, IntegrateDouble)
{
    ASSERT_TRUE(test_integrate());
}
