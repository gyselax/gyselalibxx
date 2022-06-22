#pragma once

#include <sll/spline_boundary_value.hpp>

#include "sll/view.hpp"

template <class BSplines>
class ConstantExtrapolationBoundaryValue : public SplineBoundaryValue<BSplines>
{
private:
    double m_eval_pos;

public:
    ConstantExtrapolationBoundaryValue(double eval_pos) : m_eval_pos(eval_pos) {}
    ConstantExtrapolationBoundaryValue(ConstantExtrapolationBoundaryValue const&) = delete;
    ConstantExtrapolationBoundaryValue(ConstantExtrapolationBoundaryValue&&) = delete;
    ~ConstantExtrapolationBoundaryValue() override = default;
    void operator=(ConstantExtrapolationBoundaryValue const&) = delete;
    void operator=(ConstantExtrapolationBoundaryValue&&) = delete;

    inline double operator()(double, ChunkSpan<double const, DiscreteDomain<BSplines>> spline_coef)
            const final
    {
        std::array<double, BSplines::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());
        int jmin;

        discrete_space<BSplines>().eval_basis(vals, jmin, m_eval_pos);

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines::degree() + 1; ++i) {
            y += spline_coef(DiscreteElement<BSplines>(jmin + i)) * vals(i);
        }
        return y;
    }
};
