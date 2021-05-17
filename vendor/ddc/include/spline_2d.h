#ifndef SPLINE_2D_H
#define SPLINE_2D_H
#include <memory>

#include "bsplines.h"

class Spline_2D
{
public:
    Spline_2D(const BSplines& bspl1, const BSplines& bspl2);
    ~Spline_2D();
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
            typename std::enable_if<std::is_base_of<BSplines, T1>::value>::type* = nullptr,
            typename std::enable_if<std::is_base_of<BSplines, T2>::value>::type* = nullptr>
    double eval_intern(
            double x1,
            double x2,
            const T1& bspl1,
            const T2& bspl2,
            DSpan1D& vals1,
            DSpan1D& vals2) const;
    template <
            class T1,
            typename std::enable_if<std::is_base_of<BSplines, T1>::value>::type* = nullptr,
            class T2,
            typename std::enable_if<std::is_base_of<BSplines, T2>::value>::type* = nullptr>
    void eval_array_loop(DSpan2D const& x1, DSpan2D const& x2, DSpan2D& y) const;
    std::unique_ptr<double[]> bcoef_ptr;
    DSpan2D bcoef;
    const BSplines& bspl1;
    const BSplines& bspl2;
    friend class Spline_interpolator_2D;
};

#endif // SPLINE_2D_H
