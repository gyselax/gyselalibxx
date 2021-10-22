#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/Chunk>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>

#include <sll/bsplines_uniform.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using IDimX = UniformDiscretization<DimX>;
using CoordX = Coordinate<DimX>;
using IndexX = DiscreteCoordinate<IDimX>;

class PolynomialEvaluator
{
private:
    std::vector<double> m_poly_coeff;

public:
    PolynomialEvaluator() = default;

    PolynomialEvaluator(std::vector<double> const& poly_coef) : m_poly_coeff(poly_coef) {}

    double operator()(double const x) const noexcept
    {
        return eval(x, 0);
    }

    template <class Domain>
    void operator()(ChunkSpan<double, Domain>& chunk) const
    {
        auto const& domain = chunk.domain();

        for (std::size_t i = 0; i < domain.size(); ++i) {
            chunk(i) = eval(domain.to_real(domain[i]), 0);
        }
    }

    double deriv(double const x, int const derivative) const noexcept
    {
        return eval(x, derivative);
    }

    template <class Domain>
    void deriv(ChunkSpan<double, Domain>& chunk, int const derivative) const
    {
        auto const& domain = chunk.domain();

        for (std::size_t i = 0; i < domain.size(); ++i) {
            chunk(i) = eval(domain.to_real(domain[i]), derivative);
        }
    }

private:
    double eval(double const x, int const d) const noexcept
    {
        double y = 0.0;
        for (std::size_t i = std::max(d, 0); i < m_poly_coeff.size(); ++i) {
            y += falling_factorial(i, d) * std::pow(x, (i - d)) * m_poly_coeff[i];
        }
        return y;
    }

    constexpr int falling_factorial(int x, int n) const
    {
        double c = 1.;
        if (n >= 0) {
            for (int k = 0; k < n; ++k) {
                c *= (x - k);
            }
        } else {
            for (int k = -1; k >= n; --k) {
                c /= (x - k);
            }
        }
        return c;
    }
};

class CosineEvaluator
{
    static inline constexpr double s_two_pi = 2. * M_PI;

private:
    double m_c0 = 1.;

    double m_c1 = 0.;

public:
    CosineEvaluator() = default;

    CosineEvaluator(double c0, double c1) : m_c0(c0), m_c1(c1) {}

    double operator()(double const x) const noexcept
    {
        return eval(x, 0);
    }

    template <class Domain>
    void operator()(ChunkSpan<double, Domain>& chunk) const
    {
        auto const& domain = chunk.domain();

        for (auto&& icoord : domain) {
            chunk(icoord) = eval(domain.to_real(icoord), 0);
        }
    }

    double deriv(double const x, int const derivative) const noexcept
    {
        return eval(x, derivative);
    }

    template <class Domain>
    void deriv(ChunkSpan<double, Domain>& chunk, int const derivative) const
    {
        auto const& domain = chunk.domain();

        for (auto&& icoord : domain) {
            chunk(icoord) = eval(domain.to_real(icoord), 0);
        }
    }

private:
    double eval(double const x, int const derivative) const noexcept
    {
        return std::pow(s_two_pi * m_c0, derivative)
               * std::cos(M_PI_2 * derivative + s_two_pi * (m_c0 * x + m_c1));
    }
};

TEST(SplineBuilder, Constructor)
{
    using BSplinesX = UniformBSplines<DimX, 2>;

    auto&& bsplines = BSplinesX(CoordX(0.), CoordX(0.02), 100);

    SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC> spline_builder(bsplines);
    spline_builder.interpolation_domain();
}

TEST(SplineBuilder, BuildSpline)
{
    BoundCond constexpr left_bc = DimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::HERMITE;
    BoundCond constexpr right_bc = DimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::HERMITE;
    int constexpr degree = 10;
    //     using NonUniformMeshX = NonUniformDiscretization<DimX>;
    //     using UniformMeshX = UniformDiscretization<DimX>;
    using BSplinesX = UniformBSplines<DimX, degree>;
    using SpanSplineX2 = Chunk<double, DiscreteDomain<BSplinesX>>;
    //     using NonUniformDomainX = DiscreteDomain<NonUniformMeshX>;
    //     using SpanNonUniformX = Chunk<double, DiscreteDomain<NonUniformMeshX>>;
    using SpanUniformX = Chunk<double, DiscreteDomain<IDimX>>;

    CoordX constexpr x0(0.);
    CoordX constexpr xN(1.);
    std::size_t constexpr ncells = 100;
    IndexX constexpr npoints(ncells + 1);

    // 1. Create BSplines
    BSplinesX bsplines(x0, xN, npoints);
    DiscreteDomain<BSplinesX> const
            dom_bsplines_x(bsplines, DiscreteVector<BSplinesX>(bsplines.size()));

    // 2. Create a Spline represented by a chunk over BSplines
    // The chunk is filled with garbage data, we need to initialize it
    SpanSplineX2 coef(dom_bsplines_x);

    // 3. Create a SplineBuilder over BSplines using some boundary conditions
    SplineBuilder<BSplinesX, left_bc, right_bc> spline_builder(bsplines);
    auto const& interpolation_domain = spline_builder.interpolation_domain();

    // 4. Allocate and fill a chunk over the interpolation domain
    SpanUniformX yvals(interpolation_domain);
    CosineEvaluator cosine_evaluator;
    cosine_evaluator(yvals);

    int constexpr shift = degree % 2; // shift = 0 for even order, 1 for odd order
    std::array<double, degree / 2> Sderiv_lhs_data;
    DSpan1D Sderiv_lhs(Sderiv_lhs_data.data(), Sderiv_lhs_data.size());
    for (int ii = 0; ii < Sderiv_lhs.extent(0); ++ii) {
        Sderiv_lhs(ii) = cosine_evaluator.deriv(x0, ii + shift);
    }
    std::array<double, degree / 2> Sderiv_rhs_data;
    DSpan1D Sderiv_rhs(Sderiv_rhs_data.data(), Sderiv_rhs_data.size());
    for (int ii = 0; ii < Sderiv_rhs.extent(0); ++ii) {
        Sderiv_rhs(ii) = cosine_evaluator.deriv(xN, ii + shift);
    }
    DSpan1D* deriv_l(left_bc == BoundCond::HERMITE ? &Sderiv_lhs : nullptr);
    DSpan1D* deriv_r(right_bc == BoundCond::HERMITE ? &Sderiv_rhs : nullptr);

    // 5. Finally build the spline by filling `coef`
    spline_builder(coef, yvals, deriv_l, deriv_r);

    // 6. Create a SplineEvaluator to evaluate the spline at any point in the domain of the BSplines
    SplineEvaluator spline_evaluator(bsplines, NullBoundaryValue::value, NullBoundaryValue::value);

    SpanUniformX coords_eval(interpolation_domain);
    for (auto i : interpolation_domain) {
        coords_eval(i) = interpolation_domain.to_real(i);
    }

    SpanUniformX spline_eval(interpolation_domain);
    spline_evaluator(spline_eval.span_view(), coords_eval.span_cview(), coef.span_cview());

    SpanUniformX spline_eval_deriv(interpolation_domain);
    spline_evaluator
            .deriv(spline_eval_deriv.span_view(), coords_eval.span_cview(), coef.span_cview());

    // 7. Checking errors
    double max_norm_error = 0.;
    double max_norm_error_diff = 0.;
    for (auto i : interpolation_domain) {
        auto&& x = interpolation_domain.to_real(i);

        // Compute error
        double const error = spline_eval(i) - yvals(i);
        max_norm_error = std::fmax(max_norm_error, std::fabs(error));

        // Compute error
        double const error_deriv = spline_eval_deriv(i) - cosine_evaluator.deriv(x, 1);
        max_norm_error_diff = std::fmax(max_norm_error_diff, std::fabs(error_deriv));
    }
    EXPECT_LE(max_norm_error, 1.0e-12);
}
