#ifndef SPLINE_1D_H
#define SPLINE_1D_H
#include <memory>

#include "boundary_value.hpp"
#include "bsplines.h"
#include "null_boundary_value.hpp"

class Spline_1D {
public:
    Spline_1D(const BSplines& bspl,
        const BoundaryValue& left_bc = NullBoundaryValue(),
        const BoundaryValue& right_bc = NullBoundaryValue());
    ~Spline_1D() = default;
    bool belongs_to_space(const BSplines& bspline) const;
    double eval(double x) const;
    double eval_deriv(double x) const;
    void eval_array(mdspan_1d const& x, mdspan_1d& y) const;
    void eval_array_deriv(mdspan_1d const& x, mdspan_1d& y) const;
    double integrate() const;

private:
    std::unique_ptr<double[]> bcoef_ptr;
    mdspan_1d bcoef;
    const BSplines& bspl;
    const BoundaryValue& left_bc;
    const BoundaryValue& right_bc;

    // Internal templated functions
    template <class T,
        typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
    double eval_intern_no_bcs(double x, const T& bspl, mdspan_1d& vals) const;
    template <class T,
        bool periodic,
        typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
    double eval_intern(double x, const T& bspl, mdspan_1d& vals) const;
    template <class T,
        typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
    double eval_deriv_intern(double x, const T& bspl, mdspan_1d& vals) const;
    template <class T,
        bool periodic,
        typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
    void eval_array_loop(mdspan_1d const& x, mdspan_1d& y) const;
    template <class T,
        typename std::enable_if<std::is_base_of<BSplines, T>::value>::type* = nullptr>
    void eval_array_deriv_loop(mdspan_1d const& x, mdspan_1d& y) const;

    // Friends
    friend class Spline_interpolator_1D;
    friend class Spline_interpolator_2D;
    friend int spline_1d_get_ncoeffs(Spline_1D* spl);
    friend double* spline_1d_get_coeffs(Spline_1D* spl);
};

#endif // SPLINE_1D_H
