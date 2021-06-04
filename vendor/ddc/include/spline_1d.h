#pragma once

#include <cstdint>
#include <memory>
#include <type_traits>

#include "boundary_value.h"
#include "bsplines.h"
#include "null_boundary_value.h"

class Spline1D
{
private:
    // Friends
    friend class SplineBuilder1D;
    friend class SplineBuilder2D;

private:
    std::unique_ptr<double[]> m_bcoef_ptr;
    DSpan1D m_bcoef;
    const BSplines& m_bspl;
    const BoundaryValue& m_left_bc;
    const BoundaryValue& m_right_bc;

public:
    Spline1D() = delete;
    Spline1D(
            const BSplines& bspl,
            const BoundaryValue& left_bc = NullBoundaryValue::value,
            const BoundaryValue& right_bc = NullBoundaryValue::value);
    Spline1D(const Spline1D& x) = delete;
    Spline1D(Spline1D&& x) = delete;
    ~Spline1D() = default;
    Spline1D& operator=(const Spline1D& x) = delete;
    Spline1D& operator=(Spline1D&& x) = delete;

    DSpan1D const& bcoef() const noexcept
    {
        return m_bcoef;
    }

    double const& bcoef(std::size_t i) const noexcept
    {
        return m_bcoef(i);
    }

    double& bcoef(std::size_t i) noexcept
    {
        return m_bcoef(i);
    }

    bool belongs_to_space(const BSplines& bspline) const;
    double eval(double x) const;
    double eval_deriv(double x) const;
    void eval_array(DSpan1D const x, DSpan1D y) const;
    void eval_array_deriv(DSpan1D const x, DSpan1D y) const;
    double integrate() const;

private:
    // Internal templated functions
    template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>* = nullptr>
    double eval_intern_no_bcs(double x, const T& bspl, DSpan1D& vals) const;
    template <class T, bool periodic, std::enable_if_t<std::is_base_of_v<BSplines, T>>* = nullptr>
    double eval_intern(double x, const T& bspl, DSpan1D& vals) const;
    template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>* = nullptr>
    double eval_deriv_intern(double x, const T& bspl, DSpan1D& vals) const;
    template <class T, bool periodic, std::enable_if_t<std::is_base_of_v<BSplines, T>>* = nullptr>
    void eval_array_loop(DSpan1D const& x, DSpan1D& y) const;
    template <class T, std::enable_if_t<std::is_base_of_v<BSplines, T>>* = nullptr>
    void eval_array_deriv_loop(DSpan1D const& x, DSpan1D& y) const;
};
