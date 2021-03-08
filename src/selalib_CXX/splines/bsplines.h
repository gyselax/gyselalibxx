#ifndef BSPLINES_H
#define BSPLINES_H
#include <vector>
#include <math_tools.h>

class BSplines
{
    public:
        virtual void eval_basis(double x, mdspan_1d& values, int& jmin) const = 0;
        virtual void eval_deriv(double x, mdspan_1d& derivs, int& jmin) const = 0;
        virtual void eval_basis_and_n_derivs(double x, int n, mdspan_2d& derivs, int& jmin) const = 0;
        virtual void integrate(mdspan_1d& int_vals) const = 0;
        virtual double get_knot(int idx) const = 0;
        virtual ~BSplines() {}
        static BSplines* new_bsplines(int degree, bool periodic, std::vector<double> const& breaks);
        static BSplines* new_bsplines(int degree, bool periodic, double xmin, double xmax, int ncells);
        const int degree;
        const bool periodic;
        const bool uniform;
        const bool radial;
        const int ncells;
        const int nbasis;
        const double xmin;
        const double xmax;
    protected:
        BSplines(int degree, bool periodic, bool uniform, int ncells, int nbasis, double xmin, double xmax, bool radial);
};

extern "C" {
    BSplines* get_new_bspline_uniform(int degree, bool periodic, double xmin, double xmax, int ncells);
    BSplines* get_new_bspline_non_uniform(int degree, bool periodic, double xmin, double xmax, int ncells, double* breaks_ptr, int nbreaks);
    void bsplines_free(BSplines* bspl);
    void bsplines_eval_basis(BSplines* bspl, double x, int nvals, double* vals, int* jmin);
    void bsplines_eval_deriv(BSplines* bspl, double x, int nvals, double* vals, int* jmin);
    void bsplines_eval_basis_and_n_derivs(BSplines* bspl, double x, int n, int nvals1, int nvals2, double* vals, int* jmin);
    void bsplines_integrals(BSplines* bspl, int nvals, double* int_vals);
    double bsplines_get_knot(BSplines* bspl, int idx);
}
#endif // BSPLINES_H
