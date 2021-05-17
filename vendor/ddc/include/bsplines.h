#ifndef BSPLINES_H
#define BSPLINES_H
#include <vector>

#include "view.h"

class BSplines
{
public:
    static BSplines* new_bsplines(int degree, bool periodic, std::vector<double> const& breaks);
    static BSplines* new_bsplines(int degree, bool periodic, double xmin, double xmax, int ncells);

    virtual ~BSplines() = default;

    virtual void eval_basis(double x, mdspan_1d& values, int& jmin) const = 0;
    virtual void eval_deriv(double x, mdspan_1d& derivs, int& jmin) const = 0;
    virtual void eval_basis_and_n_derivs(double x, int n, mdspan_2d& derivs, int& jmin) const = 0;
    virtual void integrals(mdspan_1d& int_vals) const = 0;
    virtual double get_knot(int idx) const = 0;

    const int degree;
    const bool uniform;
    const bool radial;
    const int nbasis;

    /// TODO: take all that from a Mesh object
    const int ncells;
    const bool periodic;
    const double xmin;
    const double xmax;
    const double length;

protected:
    BSplines(
            int degree,
            bool periodic,
            bool uniform,
            int ncells,
            int nbasis,
            double xmin,
            double xmax,
            bool radial);
};
#endif // BSPLINES_H
