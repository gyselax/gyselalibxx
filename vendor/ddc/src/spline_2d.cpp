#include <cassert>
#include <iostream>
#include <memory>

#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"
#include "spline_2d.h"

Spline2D::Spline2D(const BSplines& bspl1, const BSplines& bspl2)
    : m_bcoef_ptr(std::make_unique<double[]>(
            (bspl1.degree() + bspl1.ncells()) * (bspl2.degree() + bspl2.ncells())))
    , m_bcoef(m_bcoef_ptr.get(), bspl1.degree() + bspl1.ncells(), bspl2.degree() + bspl2.ncells())
    , m_bspl1(bspl1)
    , m_bspl2(bspl2)
{
}

bool Spline2D::belongs_to_space(const BSplines& bspline1, const BSplines& bspline2) const
{
    return (&m_bspl1 == &bspline1 && &m_bspl2 == &bspline2);
}

template <
        class T1,
        class T2,
        bool deriv1,
        bool deriv2,
        std::enable_if_t<std::is_base_of_v<BSplines, T1>>*,
        std::enable_if_t<std::is_base_of_v<BSplines, T2>>*>
double Spline2D::eval_intern(
        double x1,
        double x2,
        const T1& bspl1,
        const T2& bspl2,
        DSpan1D& vals1,
        DSpan1D& vals2) const
{
    int jmin1, jmin2;

    if constexpr (deriv1) {
        bspl1.eval_deriv(x1, vals1, jmin1);
    } else {
        bspl1.eval_basis(x1, vals1, jmin1);
    }
    if constexpr (deriv2) {
        bspl2.eval_deriv(x2, vals2, jmin2);
    } else {
        bspl2.eval_basis(x2, vals2, jmin2);
    }

    double y = 0.0;
    for (int i(0); i < bspl1.degree() + 1; ++i) {
        for (int j(0); j < bspl2.degree() + 1; ++j) {
            y += m_bcoef(jmin1 + i, jmin2 + j) * vals1(i) * vals2(j);
        }
    }
    return y;
}

double Spline2D::eval(const double x1, const double x2) const
{
    return eval_deriv<false, false>(x1, x2);
}

template <bool deriv1, bool deriv2>
double Spline2D::eval_deriv(const double x1, const double x2) const
{
    std::vector<double> values1(m_bspl1.degree() + 1);
    std::vector<double> values2(m_bspl2.degree() + 1);
    DSpan1D vals1(values1.data(), values1.size());
    DSpan1D vals2(values2.data(), values2.size());

    if (m_bspl1.is_uniform()) {
        if (m_bspl2.is_uniform()) {
            return eval_intern<UniformBSplines, UniformBSplines, deriv1, deriv2>(
                    x1,
                    x2,
                    static_cast<const UniformBSplines&>(m_bspl1),
                    static_cast<const UniformBSplines&>(m_bspl2),
                    vals1,
                    vals2);
        } else {
            return eval_intern<UniformBSplines, NonUniformBSplines, deriv1, deriv2>(
                    x1,
                    x2,
                    static_cast<const UniformBSplines&>(m_bspl1),
                    static_cast<const NonUniformBSplines&>(m_bspl2),
                    vals1,
                    vals2);
        }
    } else {
        if (m_bspl2.is_uniform()) {
            return eval_intern<NonUniformBSplines, UniformBSplines, deriv1, deriv2>(
                    x1,
                    x2,
                    static_cast<const NonUniformBSplines&>(m_bspl1),
                    static_cast<const UniformBSplines&>(m_bspl2),
                    vals1,
                    vals2);
        } else {
            return eval_intern<NonUniformBSplines, NonUniformBSplines, deriv1, deriv2>(
                    x1,
                    x2,
                    static_cast<const NonUniformBSplines&>(m_bspl1),
                    static_cast<const NonUniformBSplines&>(m_bspl2),
                    vals1,
                    vals2);
        }
    }
}

void Spline2D::eval_array(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const
{
    assert(x1.extent(0) == y.extent(0));
    assert(x1.extent(1) == y.extent(1));
    assert(x2.extent(0) == y.extent(0));
    assert(x2.extent(1) == y.extent(1));

    for (int i(0); i < x1.extent(0); ++i) {
        for (int j(0); j < x1.extent(1); ++j) {
            y(i, j) = eval(x1(i, j), x2(i, j));
        }
    }
}

template <
        class T1,
        std::enable_if_t<std::is_base_of_v<BSplines, T1>>*,
        class T2,
        std::enable_if_t<std::is_base_of_v<BSplines, T2>>*>
void Spline2D::eval_array_loop(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const
{
    std::vector<double> values1(m_bspl1.degree() + 1);
    std::vector<double> values2(m_bspl2.degree() + 1);
    DSpan1D vals1(values1.data(), values1.size());
    DSpan1D vals2(values2.data(), values2.size());

    int jmin1, jmin2;

    const T1& l_bspl1 = static_cast<const T1&>(m_bspl1);
    const T2& l_bspl2 = static_cast<const T1&>(m_bspl2);

    for (int i(0); i < x1.extent(0); ++i) {
        for (int j(0); j < x1.extent(1); ++j) {
            y(i, j) = eval_intern<
                    T1,
                    T2,
                    false,
                    false>(x1(i, j), x2(i, j), l_bspl1, l_bspl2, vals1, vals2);
        }
    }
}

template <bool deriv1, bool deriv2>
void Spline2D::eval_array_deriv(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const
{
    assert(x1.extent(0) == y.extent(0));
    assert(x1.extent(1) == y.extent(1));
    assert(x2.extent(0) == y.extent(0));
    assert(x2.extent(1) == y.extent(1));

    for (int i(0); i < x1.extent(0); ++i) {
        for (int j(0); j < x1.extent(1); ++j) {
            y(i, j) = eval_deriv<deriv1, deriv2>(x1(i, j), x2(i, j));
        }
    }
}

void Spline2D::integrate_dim(DSpan1D& y, const int dim) const
{
    assert(dim >= 0 and dim < 2);
    assert(y.extent(0) == m_bcoef.extent(1 - dim));

    const BSplines& bspline((dim == 0) ? m_bspl1 : m_bspl2);

    std::vector<double> values(m_bcoef.extent(dim));
    DSpan1D vals(values.data(), values.size());

    bspline.integrals(vals);

    if (dim == 0) {
        for (int i(0); i < y.extent(0); ++i) {
            y(i) = 0;
            for (int j(0); j < m_bcoef.extent(0); ++j) {
                y(i) += m_bcoef(j, i) * vals(j);
            }
        }
    } else {
        for (int i(0); i < y.extent(0); ++i) {
            y(i) = 0;
            for (int j(0); j < m_bcoef.extent(1); ++j) {
                y(i) += m_bcoef(i, j) * vals(j);
            }
        }
    }
}

double Spline2D::integrate() const
{
    std::vector<double> int_values(m_bcoef.extent(0));
    DSpan1D int_vals(int_values.data(), int_values.size());

    integrate_dim(int_vals, 1);

    std::vector<double> values(m_bcoef.extent(0));
    DSpan1D vals(values.data(), values.size());

    m_bspl1.integrals(vals);

    double y = 0.0;
    for (int i(0); i < m_bcoef.extent(0); ++i) {
        y += int_vals(i) * vals(i);
    }
    return y;
}
