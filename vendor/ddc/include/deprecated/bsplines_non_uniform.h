#pragma once

#include <memory>
#include <vector>

#include "bsplines.h"
#include "mdomain.h"
#include "non_uniform_mesh.h"

namespace deprecated {

class NonUniformBSplines : public BSplines
{
private:
    std::unique_ptr<double[]> m_knots;

    int m_npoints;

public:
    NonUniformBSplines() = delete;
    NonUniformBSplines(int degree, bool periodic, const std::vector<double>& breaks);
    NonUniformBSplines(const NonUniformBSplines& x) = delete;
    NonUniformBSplines(NonUniformBSplines&& x) = delete;
    virtual ~NonUniformBSplines() = default;
    NonUniformBSplines& operator=(const NonUniformBSplines& x) = delete;
    NonUniformBSplines& operator=(NonUniformBSplines&& x) = delete;

    int npoints() const noexcept
    {
        return m_npoints;
    }

    virtual void eval_basis(double x, DSpan1D& values, int& jmin) const override;
    virtual void eval_deriv(double x, DSpan1D& derivs, int& jmin) const override;
    virtual void eval_basis_and_n_derivs(double x, int n, DSpan2D& derivs, int& jmin)
            const override;
    virtual void integrals(DSpan1D& int_vals) const override;

    virtual double get_knot(int break_idx) const override
    {
        // TODO: assert break_idx >= 1 - degree
        // TODO: assert break_idx <= npoints + degree
        return m_knots[break_idx + m_degree];
    }

    bool is_uniform() const override;

private:
    int find_cell(double x) const;

    inline double& get_knot(int break_idx)
    {
        // TODO: assert break_idx >= 1 - degree
        // TODO: assert break_idx <= npoints + degree
        return m_knots[break_idx + m_degree];
    }
};

} // namespace deprecated
