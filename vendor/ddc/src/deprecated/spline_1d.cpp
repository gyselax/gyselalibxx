#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>
#include <vector>

#include "deprecated/bsplines.h"
#include "deprecated/bsplines_non_uniform.h"
#include "deprecated/bsplines_uniform.h"
#include "deprecated/spline_1d.h"

#include "boundary_value.h"

namespace deprecated {

Spline1D::Spline1D(
        const BSplines& bspl,
        const BoundaryValue& left_bc,
        const BoundaryValue& right_bc)
    : m_bcoef_ptr(std::make_unique<double[]>(bspl.degree() + bspl.ncells()))
    , m_bcoef(m_bcoef_ptr.get(), bspl.degree() + bspl.ncells())
    , m_bspl(bspl)
    , m_left_bc(left_bc)
    , m_right_bc(right_bc)
{
}

bool Spline1D::belongs_to_space(const BSplines& bspline) const
{
    return &m_bspl == &bspline;
}

template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>*>
double Spline1D::eval_intern_no_bcs(double x, const T& bspl, DSpan1D& vals) const
{
    int jmin;

    bspl.eval_basis(x, vals, jmin);

    double y = 0.0;
    for (int i(0); i < bspl.degree() + 1; ++i) {
        y += m_bcoef(jmin + i) * vals(i);
    }
    return y;
}

template <class T, bool periodic, std::enable_if_t<std::is_base_of_v<BSplines, T>>*>
double Spline1D::eval_intern(double x, const T& bspl, DSpan1D& vals) const
{
    if constexpr (periodic) {
        if (x < bspl.xmin() || x > bspl.xmax())
            [[unlikely]]
            {
                x -= std::floor((x - bspl.xmin()) / bspl.length()) * bspl.length();
            }
    } else {
        if (x < bspl.xmin())
            [[unlikely]]
            {
                return m_left_bc(x);
            }
        if (x > bspl.xmax())
            [[unlikely]]
            {
                return m_right_bc(x);
            }
    }
    return eval_intern_no_bcs<T>(x, bspl, vals);
}

template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>*>
double Spline1D::eval_deriv_intern(double x, const T& bspl, DSpan1D& vals) const
{
    int jmin;

    bspl.eval_deriv(x, vals, jmin);

    double y = 0.0;
    for (int i(0); i < bspl.degree() + 1; ++i) {
        y += m_bcoef(jmin + i) * vals(i);
    }
    return y;
}

double Spline1D::eval(double x) const
{
    std::vector<double> values(m_bspl.degree() + 1);
    DSpan1D vals(values.data(), values.size());

    if (m_bspl.is_uniform()) {
        if (m_bspl.is_periodic()) {
            return eval_intern<
                    UniformBSplines,
                    true>(x, static_cast<const UniformBSplines&>(m_bspl), vals);
        } else {
            return eval_intern<
                    UniformBSplines,
                    false>(x, static_cast<const UniformBSplines&>(m_bspl), vals);
        }
    } else {
        if (m_bspl.is_periodic()) {
            return eval_intern<
                    NonUniformBSplines,
                    true>(x, static_cast<const NonUniformBSplines&>(m_bspl), vals);
        } else {
            return eval_intern<
                    NonUniformBSplines,
                    false>(x, static_cast<const NonUniformBSplines&>(m_bspl), vals);
        }
    }
}

double Spline1D::eval_deriv(double x) const
{
    std::vector<double> values(m_bspl.degree() + 1);
    DSpan1D vals(values.data(), values.size());

    if (m_bspl.is_uniform())
        return eval_deriv_intern<
                UniformBSplines>(x, static_cast<const UniformBSplines&>(m_bspl), vals);
    else
        return eval_deriv_intern<
                NonUniformBSplines>(x, static_cast<const NonUniformBSplines&>(m_bspl), vals);
}

template <class T, bool periodic, std::enable_if_t<std::is_base_of_v<BSplines, T>>*>
void Spline1D::eval_array_loop(DSpan1D const& x, DSpan1D& y) const
{
    const T& l_bspl = static_cast<const T&>(m_bspl);

    assert(x.extent(0) == y.extent(0));
    std::vector<double> values(l_bspl.degree() + 1);
    DSpan1D vals(values.data(), values.size());

    for (int i(0); i < x.extent(0); ++i) {
        y(i) = eval_intern<T, periodic>(x(i), l_bspl, vals);
    }
}

template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>*>
void Spline1D::eval_array_deriv_loop(DSpan1D const& x, DSpan1D& y) const
{
    const T& l_bspl = static_cast<const T&>(m_bspl);

    assert(x.extent(0) == y.extent(0));
    std::vector<double> values(l_bspl.degree() + 1);
    DSpan1D vals(values.data(), values.size());

    for (int i(0); i < x.extent(0); ++i) {
        y(i) = eval_deriv_intern<T>(x(i), l_bspl, vals);
    }
}

void Spline1D::eval_array(DSpan1D const x, DSpan1D y) const
{
    if (m_bspl.is_uniform()) {
        if (m_bspl.is_periodic()) {
            return eval_array_loop<UniformBSplines, true>(x, y);
        } else {
            return eval_array_loop<UniformBSplines, false>(x, y);
        }
    } else {
        if (m_bspl.is_periodic()) {
            return eval_array_loop<NonUniformBSplines, true>(x, y);
        } else {
            return eval_array_loop<NonUniformBSplines, false>(x, y);
        }
    }
}

void Spline1D::eval_array_deriv(DSpan1D const x, DSpan1D y) const
{
    if (m_bspl.is_uniform())
        eval_array_deriv_loop<UniformBSplines>(x, y);
    else
        eval_array_deriv_loop<NonUniformBSplines>(x, y);
}

double Spline1D::integrate() const
{
    std::vector<double> values(m_bcoef.extent(0));
    DSpan1D vals(values.data(), values.size());

    m_bspl.integrals(vals);

    double y = 0.0;
    for (int i(0); i < m_bcoef.extent(0); ++i) {
        y += m_bcoef(i) * vals(i);
    }
    return y;
}

} // namespace deprecated
