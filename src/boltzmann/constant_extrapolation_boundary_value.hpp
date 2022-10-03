#pragma once

#include <sll/spline_boundary_value.hpp>

#include "sll/view.hpp"

template <class BSplines>
class ConstantExtrapolationBoundaryValue : public SplineBoundaryValue<BSplines>
{
private:
    double m_eval_pos;

public:
    explicit ConstantExtrapolationBoundaryValue(double eval_pos) : m_eval_pos(eval_pos) {}

    ~ConstantExtrapolationBoundaryValue() override = default;

    double operator()(double, ChunkSpan<double const, DiscreteDomain<BSplines>> const spline_coef)
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
