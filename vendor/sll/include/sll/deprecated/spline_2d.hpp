#pragma once

#include <iosfwd>
#include <memory>
#include <type_traits>

#include "sll/view.hpp"

namespace deprecated {
class BSplines;

class Spline2D
{
private:
    friend class SplineBuilder2D;

private:
    std::unique_ptr<double[]> m_bcoef_ptr;
    DSpan2D m_bcoef;
    const BSplines& m_bspl1;
    const BSplines& m_bspl2;

public:
    Spline2D() = delete;
    Spline2D(const BSplines& bspl1, const BSplines& bspl2);
    Spline2D(const Spline2D& x) = delete;
    Spline2D(Spline2D&& x) = delete;
    ~Spline2D() = default;
    Spline2D& operator=(const Spline2D& x) = delete;
    Spline2D& operator=(Spline2D&& x) = delete;

    DSpan2D const& bcoef() const noexcept
    {
        return m_bcoef;
    }

    double& bcoef(std::size_t i, std::size_t j) const noexcept
    {
        return m_bcoef(i, j);
    }

    bool belongs_to_space(const BSplines& bspline1, const BSplines& bspline2) const;
    double eval(const double x1, const double x2) const;
    template <bool deriv1, bool deriv2>
    double eval_deriv(const double x1, const double x2) const;
    void eval_array(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const;
    template <bool deriv1, bool deriv2>
    void eval_array_deriv(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const;
    void integrate_dim(DSpan1D& y, const int dim) const;
    double integrate() const;

private:
    template <
            class T1,
            class T2,
            bool deriv1,
            bool deriv2,
            std::enable_if_t<std::is_base_of_v<BSplines, T1>>* = nullptr,
            std::enable_if_t<std::is_base_of_v<BSplines, T2>>* = nullptr>
    double eval_intern(
            double x1,
            double x2,
            const T1& bspl1,
            const T2& bspl2,
            DSpan1D& vals1,
            DSpan1D& vals2) const;
    template <
            class T1,
            std::enable_if_t<std::is_base_of_v<BSplines, T1>>* = nullptr,
            class T2,
            std::enable_if_t<std::is_base_of_v<BSplines, T2>>* = nullptr>
    void eval_array_loop(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const;
};

} // namespace deprecated
