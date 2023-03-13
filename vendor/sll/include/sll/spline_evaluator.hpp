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

    using tag_type = typename BSplinesType::tag_type;

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
            ddc::Coordinate<tag_type> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals = as_span(values);

        return eval(coord_eval, spline_coef, vals);
    }

    template <class Domain>
    void operator()(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<const ddc::Coordinate<tag_type>, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals = as_span(values);

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval(coords_eval(i), spline_coef, vals);
        }
    }

    double deriv(
            ddc::Coordinate<tag_type> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals = as_span(values);

        return eval_no_bc(coord_eval, spline_coef, vals, eval_deriv_type());
    }

    template <class Domain>
    void deriv(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<const ddc::Coordinate<tag_type>, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        std::array<double, bsplines_type::degree() + 1> values;
        DSpan1D const vals = as_span(values);

        for (auto i : coords_eval.domain()) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, vals, eval_deriv_type());
        }
    }

    double integrate(
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef) const
    {
        Chunk<double, ddc::DiscreteDomain<BSplinesType>> values(spline_coef.domain());

        ddc::discrete_space<bsplines_type>().integrals(values);

        return ddc::transform_reduce(
                ddc::policies::parallel_host,
                spline_coef.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                [&](ddc::DiscreteElement<BSplinesType> const ibspl) {
                    return spline_coef(ibspl) * values(ibspl);
                });
    }

private:
    double eval(
            ddc::Coordinate<tag_type> coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef,
            DSpan1D const vals) const
    {
        if constexpr (bsplines_type::is_periodic()) {
            if (coord_eval < ddc::discrete_space<bsplines_type>().rmin()
                || coord_eval > ddc::discrete_space<bsplines_type>().rmax()) {
                coord_eval -= std::floor(
                                      (coord_eval - ddc::discrete_space<bsplines_type>().rmin())
                                      / ddc::discrete_space<bsplines_type>().length())
                              * ddc::discrete_space<bsplines_type>().length();
            }
        } else {
            if (coord_eval < ddc::discrete_space<bsplines_type>().rmin()) {
                return m_left_bc(coord_eval, spline_coef);
            }
            if (coord_eval > ddc::discrete_space<bsplines_type>().rmax()) {
                return m_right_bc(coord_eval, spline_coef);
            }
        }
        return eval_no_bc(coord_eval, spline_coef, vals, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            ddc::Coordinate<tag_type> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType>> const spline_coef,
            DSpan1D const vals,
            EvalType const) const
    {
        static_assert(
                std::is_same_v<EvalType, eval_type> || std::is_same_v<EvalType, eval_deriv_type>);
        DiscreteElement<BSplinesType> jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = ddc::discrete_space<bsplines_type>().eval_basis(vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_type>) {
            jmin = ddc::discrete_space<bsplines_type>().eval_deriv(vals, coord_eval);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            y += spline_coef(jmin + i) * vals(i);
        }
        return y;
    }
};
