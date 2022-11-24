#pragma once

#include <array>
#include <cmath>

#include <ddc/ddc.hpp>

#include "sll/spline_boundary_value.hpp"
#include "sll/view.hpp"

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
    SplineBoundaryValue<BSplinesType> const& m_left_bc;

    SplineBoundaryValue<BSplinesType> const& m_right_bc;

public:
    SplineEvaluator() = delete;

    explicit SplineEvaluator(
            SplineBoundaryValue<BSplinesType> const& left_bc,
            SplineBoundaryValue<BSplinesType> const& right_bc)
        : m_left_bc(left_bc)
        , m_right_bc(right_bc)
    {
    }

    SplineEvaluator(SplineEvaluator const& x) = default;

    SplineEvaluator(SplineEvaluator&& x) = default;

    ~SplineEvaluator() = default;

    SplineEvaluator& operator=(SplineEvaluator const& x) = default;

    SplineEvaluator& operator=(SplineEvaluator&& x) = default;

    double operator()(
            double const coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());

        return eval(coord_eval, spline_coef, vals);
    }

    template <class Domain>
    void operator()(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<double const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval(coords_eval(i), spline_coef, vals);
        }
    }

    double deriv(
            double const coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());

        return eval_no_bc(coord_eval, spline_coef, vals, eval_deriv_type());
    }

    template <class Domain>
    void deriv(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<double const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals(values.data(), values.size());

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, vals, eval_deriv_type());
        }
    }

    double integrate(ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        Chunk<double, DiscreteDomain<BSplinesType>> values(spline_coef.domain());

        discrete_space<bsplines_type>().integrals(values);

        return transform_reduce(
                policies::parallel_host,
                spline_coef.domain(),
                0.0,
                reducer::sum<double>(),
                [&](DiscreteElement<BSplinesType> const ibspl) {
                    return spline_coef(ibspl) * values(ibspl);
                });
    }

private:
    double eval(
            double coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef,
            DSpan1D const vals) const
    {
        if constexpr (bsplines_type::is_periodic()) {
            if (coord_eval < discrete_space<bsplines_type>().rmin()
                || coord_eval > discrete_space<bsplines_type>().rmax()) {
                coord_eval -= std::floor(
                                      (coord_eval - discrete_space<bsplines_type>().rmin())
                                      / discrete_space<bsplines_type>().length())
                              * discrete_space<bsplines_type>().length();
            }
        } else {
            if (coord_eval < discrete_space<bsplines_type>().rmin()) {
                return m_left_bc(coord_eval, spline_coef);
            }
            if (coord_eval > discrete_space<bsplines_type>().rmax()) {
                return m_right_bc(coord_eval, spline_coef);
            }
        }
        return eval_no_bc(coord_eval, spline_coef, vals, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            double const coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType>> const spline_coef,
            DSpan1D const vals,
            EvalType const) const
    {
        static_assert(
                std::is_same_v<EvalType, eval_type> || std::is_same_v<EvalType, eval_deriv_type>);
        DiscreteElement<BSplinesType> jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = discrete_space<bsplines_type>().eval_basis(vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_type>) {
            jmin = discrete_space<bsplines_type>().eval_deriv(vals, coord_eval);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            y += spline_coef(jmin + i) * vals(i);
        }
        return y;
    }
};
