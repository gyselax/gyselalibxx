#ifndef SPLINE_2D_H
#define SPLINE_2D_H
#include "bsplines.h"

class Spline_2D
{
    public:
        Spline_2D(BSplines* bspl1, BSplines* bspl2);
        ~Spline_2D();
        bool belongs_to_space(const BSplines* bspline1, const BSplines* bspline2) const;
        double eval(const double x1, const double x2) const;
        double eval_deriv(const double x1, const double x2, const bool deriv1, const bool deriv2) const;
        void eval_array(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y) const;
        void eval_array_deriv(mdspan_2d const& x1, mdspan_2d const& x2, mdspan_2d& y, const bool deriv1, const bool deriv2) const;
        void integrate_dim(mdspan_1d& y, const int dim) const;
        double integrate() const;
    private:
        double* bcoef_ptr;
        mdspan_2d bcoef;
        const BSplines* bspl1;
        const BSplines* bspl2;
        friend class Spline_interpolator_2D;
};

extern "C"
{
    Spline_2D* spline_2d_new(BSplines* bspl1, BSplines* bspl2);
    void spline_2d_free(Spline_2D* spl);
    double spline_2d_eval(Spline_2D* spl, double* x1, double* x2);
    double spline_2d_eval_deriv(Spline_2D* spl, double* x1, double* x2,
            bool deriv_x1, bool deriv_x2);
    void spline_2d_eval_array(Spline_2D* spl,
            double* x1_ptr, int nx1_1, int nx1_2,
            double* x2_ptr, int nx2_1, int nx2_2,
            double* y_ptr,  int ny1,   int ny2);
    void spline_2d_eval_array_deriv(Spline_2D* spl,
            double* x1_ptr, int nx1_1, int nx1_2,
            double* x2_ptr, int nx2_1, int nx2_2,
            bool deriv_x1,  bool deriv_x2,
            double* y_ptr,  int ny1, int ny2);
    void spline_2d_integrate_dim(Spline_2D* spl, double* y_ptr, int ny, const int dim);
    double spline_2d_integrate(Spline_2D* spl);
}

#endif // SPLINE_2D_H
