#pragma once
#include <sll/spline_boundary_value.hpp>

template <class BSplines>
class ConstantExtrapolationBoundaryValue : public SplineBoundaryValue<BSplines>
{
    ConstantExtrapolationBoundaryValue(ConstantExtrapolationBoundaryValue const&) = delete;
    ConstantExtrapolationBoundaryValue(ConstantExtrapolationBoundaryValue&&) = delete;
    void operator=(ConstantExtrapolationBoundaryValue const&) = delete;
    void operator=(ConstantExtrapolationBoundaryValue&&) = delete;

private:
    double m_eval_pos;

public:
    ConstantExtrapolationBoundaryValue(double eval_pos) : m_eval_pos(eval_pos) {}
    ~ConstantExtrapolationBoundaryValue() override = default;

    inline double operator()(
            double coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplines>> spline_coef) const final
    {
        std::array<double, BSplines::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());
        int jmin;

        discretization<BSplines>().eval_basis(vals, jmin, coord_eval);

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines::degree() + 1; ++i) {
            y += spline_coef(DiscreteCoordinate<BSplines>(jmin + i)) * vals(i);
        }
        return y;
    }
};
