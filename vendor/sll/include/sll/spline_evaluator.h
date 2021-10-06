#pragma once

#include <array>
#include <cmath>

#include <ddc/BlockSpan>

#include "sll/boundary_value.h"

template <class BSplinesType>
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
    using bsplines_type = BSplinesType;

private:
    bsplines_type const& m_bsplines;

    BoundaryValue const& m_left_bc;

    BoundaryValue const& m_right_bc;

public:
    SplineEvaluator() = delete;

    explicit SplineEvaluator(
            BSplinesType const& bsplines,
            BoundaryValue const& left_bc,
            BoundaryValue const& right_bc)
        : m_bsplines(bsplines)
        , m_left_bc(left_bc)
        , m_right_bc(right_bc)
    {
    }

    SplineEvaluator(const SplineEvaluator& x) = default;

    SplineEvaluator(SplineEvaluator&& x) = default;

    ~SplineEvaluator() = default;

    SplineEvaluator& operator=(const SplineEvaluator& x) = default;

    SplineEvaluator& operator=(SplineEvaluator&& x) = default;

    double operator()(
            double coord_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        return eval(coord_eval, spline_coef, vals);
    }

    template <class Domain>
    void operator()(
            BlockSpan<double, Domain> const& spline_eval,
            BlockSpan<double const, Domain> const& coords_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval(coords_eval(i), spline_coef, vals);
        }
    }

    double deriv(
            double coord_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        return eval_no_bc(coord_eval, spline_coef, vals, eval_deriv_type());
    }

    template <class Domain>
    void deriv(
            BlockSpan<double, Domain> const& spline_eval,
            BlockSpan<double const, Domain> const& coords_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                vals(values.data(), values.size());

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, vals, eval_deriv_type());
        }
    }

    double integrate(BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef) const
    {
        Block<double, ProductMDomain<BSplinesType>> values(spline_coef.domain());
        auto vals = values.allocation_mdspan();

        m_bsplines.integrals(vals);

        double y = 0.;
        for (MCoord<BSplinesType> ibs : spline_coef.domain()) {
            y += spline_coef(ibs) * values(ibs);
        }
        return y;
    }

private:
    double eval(
            double coord_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& vals) const
    {
        if constexpr (bsplines_type::is_periodic()) {
            if (coord_eval < m_bsplines.rmin() || coord_eval > m_bsplines.rmax()) {
                coord_eval -= std::floor((coord_eval - m_bsplines.rmin()) / m_bsplines.length())
                              * m_bsplines.length();
            }
        } else {
            if (coord_eval < m_bsplines.rmin()) {
                return m_left_bc(coord_eval);
            }
            if (coord_eval > m_bsplines.rmax()) {
                return m_right_bc(coord_eval);
            }
        }
        return eval_no_bc(coord_eval, spline_coef, vals, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            double coord_eval,
            BlockSpan<double const, ProductMDomain<BSplinesType>> const& spline_coef,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& vals,
            EvalType) const
    {
        static_assert(
                std::is_same_v<EvalType, eval_type> || std::is_same_v<EvalType, eval_deriv_type>);
        int jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            m_bsplines.eval_basis(coord_eval, vals, jmin);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_type>) {
            m_bsplines.eval_deriv(coord_eval, vals, jmin);
        }

        double y = 0.0;
        for (int i(0); i < bsplines_type::degree() + 1; ++i) {
            y += spline_coef(jmin + i) * vals(i);
        }
        return y;
    }
};
