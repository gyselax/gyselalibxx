#pragma once

#include <array>
#include <cmath>

#include <ddc/BlockSpan>

#include "sll/boundary_value.h"
#include "sll/bsplines.h"

template <class BlockSplineType>
class SplineEvaluator
{
private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

public:
    using spline_type = BlockSplineType;

    using bsplines_type = typename spline_type::bsplines_type;

private:
    spline_type const& m_spline;

    BoundaryValue const& m_left_bc;

    BoundaryValue const& m_right_bc;

public:
    SplineEvaluator() = delete;

    explicit SplineEvaluator(
            BlockSplineType const& spline,
            BoundaryValue const& left_bc,
            BoundaryValue const& right_bc)
        : m_spline(spline)
        , m_left_bc(left_bc)
        , m_right_bc(right_bc)
    {
    }

    SplineEvaluator(const SplineEvaluator& x) = default;

    SplineEvaluator(SplineEvaluator&& x) = default;

    ~SplineEvaluator() = default;

    SplineEvaluator& operator=(const SplineEvaluator& x) = default;

    SplineEvaluator& operator=(SplineEvaluator&& x) = default;

    double operator()(double x) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        return eval(x, vals);
    }

    template <class Domain>
    void operator()(BlockSpan<Domain, double>& block_mesh) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        auto const& domain = block_mesh.domain();

        for (std::size_t i = 0; i < domain.size(); ++i) {
            block_mesh(i) = eval(domain.to_real(i), vals);
        }
    }

    double deriv(double x) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        return eval_no_bc(x, vals, eval_deriv_type());
    }

    template <class Domain>
    void deriv(BlockSpan<Domain, double>& block_mesh) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        auto const& domain = block_mesh.domain();

        for (std::size_t i = 0; i < domain.size(); ++i) {
            block_mesh(i) = eval_no_bc(domain.to_real(i), vals, eval_deriv_type());
        }
    }

private:
    double eval(double x, std::experimental::mdspan<double, std::experimental::dextents<1>>& vals)
            const
    {
        if constexpr (bsplines_type::is_periodic()) {
            if (x < m_spline.bsplines().rmin() || x > m_spline.bsplines().rmax()) {
                x -= std::floor((x - m_spline.bsplines().rmin()) / m_spline.bsplines().length())
                     * m_spline.bsplines().length();
            }
        } else {
            if (x < m_spline.bsplines().rmin()) {
                return m_left_bc(x);
            }
            if (x > m_spline.bsplines().rmax()) {
                return m_right_bc(x);
            }
        }
        return eval_no_bc(x, vals, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            double x,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& vals,
            EvalType) const
    {
        static_assert(
                std::is_same_v<EvalType, eval_type> || std::is_same_v<EvalType, eval_deriv_type>);
        int jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            m_spline.bsplines().eval_basis(x, vals, jmin);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_type>) {
            m_spline.bsplines().eval_deriv(x, vals, jmin);
        }

        double y = 0.0;
        for (int i(0); i < m_spline.bsplines().degree() + 1; ++i) {
            y += m_spline(jmin + i) * vals(i);
        }
        return y;
    }
};
