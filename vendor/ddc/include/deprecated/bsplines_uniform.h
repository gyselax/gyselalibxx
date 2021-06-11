#pragma once

#include "bsplines.h"
#include "mdomain.h"

namespace deprecated {

class UniformBSplines : public BSplines
{
private:
    double m_inv_dx;

    double m_dx;

public:
    UniformBSplines() = delete;
    template <class Tag>
    UniformBSplines(int degree, const MDomainImpl<UniformMesh<Tag>>& dom)
        : UniformBSplines(
                degree,
                Tag::PERIODIC,
                dom.rmin(),
                dom.rmax() - dom.mesh().step(),
                dom.size() - 1)
    {
    }
    UniformBSplines(int degree, bool periodic, double xmin, double xmax, int ncells);
    UniformBSplines(const UniformBSplines& x) = delete;
    UniformBSplines(UniformBSplines&& x) = delete;
    virtual ~UniformBSplines() = default;
    UniformBSplines& operator=(const UniformBSplines& x) = delete;
    UniformBSplines& operator=(UniformBSplines&& x) = delete;
    virtual inline void eval_basis(double x, DSpan1D& values, int& jmin) const override
    {
        return eval_basis(x, values, jmin, m_degree);
    }
    virtual void eval_deriv(double x, DSpan1D& derivs, int& jmin) const override;
    virtual void eval_basis_and_n_derivs(double x, int n, DSpan2D& derivs, int& jmin)
            const override;
    virtual void integrals(DSpan1D& int_vals) const override;

    virtual double get_knot(int idx) const override
    {
        return m_xmin + idx * m_dx;
    }

    bool is_uniform() const override;

private:
    void eval_basis(double x, DSpan1D& values, int& jmin, int degree) const;
    void get_icell_and_offset(double x, int& icell, double& offset) const;
};

} // namespace deprecated
