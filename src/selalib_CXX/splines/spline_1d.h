#ifndef SPLINE_1D_H
#define SPLINE_1D_H
#include "bsplines.h"

class Spline_1D;

extern "C"
{
    Spline_1D* spline_1d_new(BSplines* bspl);
    void spline_1d_free(Spline_1D* spl);
    double spline_1d_eval(Spline_1D* spl, double x);
    double spline_1d_eval_deriv(Spline_1D* spl, double x);
    void spline_1d_eval_array(Spline_1D* spl, double* x, int nx, double* y, int ny);
    void spline_1d_eval_array_deriv(Spline_1D* spl, double* x, int nx, double* y, int ny);
    double spline_1d_integrate(Spline_1D* spl);
    int spline_1d_get_ncoeffs(Spline_1D* spl);
    double* spline_1d_get_coeffs(Spline_1D* spl);
}

class Spline_1D
{
    public:
        Spline_1D(const BSplines* bspl);
        ~Spline_1D();
        bool belongs_to_space(const BSplines* bspline) const;
        double eval(double x) const;
        double eval_deriv(double x) const;
        void eval_array(mdspan_1d const& x, mdspan_1d& y) const;
        void eval_array_deriv(mdspan_1d const& x, mdspan_1d& y) const;
        double integrate() const;
    private:
        double* bcoef_ptr;
        mdspan_1d bcoef;
        const BSplines* bspl;
        friend class Spline_interpolator_1D;
        friend class Spline_interpolator_2D;
        friend int spline_1d_get_ncoeffs(Spline_1D* spl);
        friend double* spline_1d_get_coeffs(Spline_1D* spl);
};

#endif // SPLINE_1D_H
