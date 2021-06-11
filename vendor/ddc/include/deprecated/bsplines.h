#pragma once

#include "view.h"

namespace deprecated {

class BSplines
{
protected:
    int m_degree;
    bool m_radial;
    int m_nbasis;

    /// TODO: take all that from a Mesh object
    int m_ncells;
    bool m_periodic;
    double m_xmin;
    double m_xmax;
    double m_length;

public:
    virtual ~BSplines() = default;

    constexpr inline int degree() const
    {
        return m_degree;
    }

    constexpr inline bool radial() const
    {
        return m_radial;
    }

    constexpr inline int nbasis() const
    {
        return m_nbasis;
    }

    constexpr inline int ncells() const
    {
        return m_ncells;
    }

    constexpr inline bool is_periodic() const
    {
        return m_periodic;
    }

    constexpr inline double xmin() const
    {
        return m_xmin;
    }

    constexpr inline double xmax() const
    {
        return m_xmax;
    }

    constexpr inline double length() const
    {
        return m_length;
    }

    virtual void eval_basis(double x, DSpan1D& values, int& jmin) const = 0;

    virtual void eval_deriv(double x, DSpan1D& derivs, int& jmin) const = 0;

    virtual void eval_basis_and_n_derivs(double x, int n, DSpan2D& derivs, int& jmin) const = 0;

    virtual void integrals(DSpan1D& int_vals) const = 0;

    virtual double get_knot(int idx) const = 0;

    virtual bool is_uniform() const = 0;

protected:
    BSplines(
            int degree,
            bool periodic,
            int ncells,
            int nbasis,
            double xmin,
            double xmax,
            bool radial);
};

} // namespace deprecated
