#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include <sll/spline_boundary_value.hpp>
#include <sll/view.hpp>

template <class BSplinesType1, class BSplinesType2>
class SplineEvaluator2D
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
    using bsplines_type1 = BSplinesType1;
    using bsplines_type2 = BSplinesType2;
    using Dim1 = typename BSplinesType1::tag_type;
    using Dim2 = typename BSplinesType2::tag_type;

private:
    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_left_bc_1;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_right_bc_1;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_left_bc_2;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_right_bc_2;

public:
    SplineEvaluator2D() = delete;

    explicit SplineEvaluator2D(
            SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& left_bc_1,
            SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& right_bc_1,
            SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& left_bc_2,
            SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& right_bc_2)
        : m_left_bc_1(left_bc_1)
        , m_right_bc_1(right_bc_1)
        , m_left_bc_2(left_bc_2)
        , m_right_bc_2(right_bc_2)
    {
    }

    SplineEvaluator2D(SplineEvaluator2D const& x) = default;

    SplineEvaluator2D(SplineEvaluator2D&& x) = default;

    ~SplineEvaluator2D() = default;

    SplineEvaluator2D& operator=(SplineEvaluator2D const& x) = default;

    SplineEvaluator2D& operator=(SplineEvaluator2D&& x) = default;

    double operator()(
            Coordinate<Dim1, Dim2> const& coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        return eval(coord_eval, spline_coef, vals1, vals2);
    }

    template <class Domain>
    void operator()(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2(values2.data(), values2.size());

        for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval(coords_eval(i), spline_coef, vals1, vals2);
        });
    }

    double deriv_dim_1(
            Coordinate<Dim1, Dim2> const& coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        return eval_no_bc(
                get<Dim1>(coord_eval),
                get<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_deriv_type(),
                eval_type());
    }

    double deriv_dim_2(
            Coordinate<Dim1, Dim2> const& coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        return eval_no_bc(
                get<Dim1>(coord_eval),
                get<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_type(),
                eval_deriv_type());
    }

    double deriv_1_and_2(
            Coordinate<Dim1, Dim2> const& coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        return eval_no_bc(
                get<Dim1>(coord_eval),
                get<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_deriv_type(),
                eval_deriv_type());
    }

    template <class Domain>
    void deriv_dim_1(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    get<Dim1>(coords_eval(i)),
                    get<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_deriv_type(),
                    eval_type());
        });
    }

    template <class Domain>
    void deriv_dim_2(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    get<Dim1>(coords_eval(i)),
                    get<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_type(),
                    eval_deriv_type());
        });
    }

    template <class Domain>
    void deriv_dim_1_and_2(
            ChunkSpan<double, Domain> const spline_eval,
            ChunkSpan<Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef)
            const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D vals1(values1.data(), values1.size());
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D vals2(values2.data(), values2.size());

        for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    get<Dim1>(coords_eval(i)),
                    get<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_deriv_type(),
                    eval_deriv_type());
        });
    }

    double integrate(ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const
                             spline_coef) const
    {
        Chunk<double, DiscreteDomain<BSplinesType1>> values1(spline_coef.domain());
        DSpan1D vals1 = values1.allocation_mdspan();
        Chunk<double, DiscreteDomain<BSplinesType2>> values2(spline_coef.domain());
        DSpan1D vals2 = values2.allocation_mdspan();

        discrete_space<bsplines_type1>().integrals(vals1);
        discrete_space<bsplines_type2>().integrals(vals2);

        return transform_reduce(
                policies::parallel_host,
                spline_coef.domain(),
                0.0,
                reducer::sum<double>(),
                [&](DiscreteElement<BSplinesType1, BSplinesType2> const i) {
                    return spline_coef(i) * values1(select<BSplinesType1>(i))
                           * values2(select<BSplinesType2>(i));
                });
    }

private:
    double eval(
            Coordinate<Dim1, Dim2> const& coord_eval,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef,
            DSpan1D const vals1,
            DSpan1D const vals2) const
    {
        double coord_eval1 = get<Dim1>(coord_eval);
        double coord_eval2 = get<Dim2>(coord_eval);
        if constexpr (bsplines_type1::is_periodic()) {
            if (coord_eval1 < discrete_space<bsplines_type1>().rmin()
                || coord_eval1 > discrete_space<bsplines_type1>().rmax()) {
                coord_eval1 -= std::floor(
                                       (coord_eval1 - discrete_space<bsplines_type1>().rmin())
                                       / discrete_space<bsplines_type1>().length())
                               * discrete_space<bsplines_type1>().length();
            }
        } else {
            if (coord_eval1 < discrete_space<bsplines_type1>().rmin()) {
                return m_left_bc_1(coord_eval1, coord_eval2, spline_coef);
            }
            if (coord_eval1 > discrete_space<bsplines_type1>().rmax()) {
                return m_right_bc_1(coord_eval1, coord_eval2, spline_coef);
            }
        }
        if constexpr (bsplines_type2::is_periodic()) {
            if (coord_eval2 < discrete_space<bsplines_type2>().rmin()
                || coord_eval2 > discrete_space<bsplines_type2>().rmax()) {
                coord_eval2 -= std::floor(
                                       (coord_eval2 - discrete_space<bsplines_type2>().rmin())
                                       / discrete_space<bsplines_type2>().length())
                               * discrete_space<bsplines_type2>().length();
            }
        } else {
            if (coord_eval2 < discrete_space<bsplines_type2>().rmin()) {
                return m_left_bc_2(coord_eval1, coord_eval2, spline_coef);
            }
            if (coord_eval2 > discrete_space<bsplines_type2>().rmax()) {
                return m_right_bc_2(coord_eval1, coord_eval2, spline_coef);
            }
        }
        return eval_no_bc(
                coord_eval1,
                coord_eval2,
                spline_coef,
                vals1,
                vals2,
                eval_type(),
                eval_type());
    }

    template <class EvalType1, class EvalType2>
    double eval_no_bc(
            double const coord_eval1,
            double const coord_eval2,
            ChunkSpan<double const, DiscreteDomain<BSplinesType1, BSplinesType2>> const spline_coef,
            DSpan1D const vals1,
            DSpan1D const vals2,
            EvalType1 const,
            EvalType2 const) const
    {
        static_assert(
                std::is_same_v<EvalType1, eval_type> || std::is_same_v<EvalType1, eval_deriv_type>);
        static_assert(
                std::is_same_v<EvalType2, eval_type> || std::is_same_v<EvalType2, eval_deriv_type>);
        DiscreteElement<BSplinesType1> jmin1;
        DiscreteElement<BSplinesType2> jmin2;

        if constexpr (std::is_same_v<EvalType1, eval_type>) {
            jmin1 = discrete_space<bsplines_type1>().eval_basis(vals1, coord_eval1);
        } else if constexpr (std::is_same_v<EvalType1, eval_deriv_type>) {
            jmin1 = discrete_space<bsplines_type1>().eval_deriv(vals1, coord_eval1);
        }
        if constexpr (std::is_same_v<EvalType2, eval_type>) {
            jmin2 = discrete_space<bsplines_type2>().eval_basis(vals2, coord_eval2);
        } else if constexpr (std::is_same_v<EvalType2, eval_deriv_type>) {
            jmin2 = discrete_space<bsplines_type2>().eval_deriv(vals2, coord_eval2);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < bsplines_type2::degree() + 1; ++j) {
                y += spline_coef(jmin1 + i, jmin2 + j) * vals1(i) * vals2(j);
            }
        }
        return y;
    }
};
