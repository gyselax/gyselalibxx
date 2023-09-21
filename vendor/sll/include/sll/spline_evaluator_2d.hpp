#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include <sll/spline_boundary_value.hpp>
#include <sll/view.hpp>

/**
 * @brief Define an evaluator 2D on B-splines.
 */
template <class BSplinesType1, class BSplinesType2>
class SplineEvaluator2D
{
private:
    /**
     * @brief Tag to indicate that the value of the spline should be evaluated.
     */
    struct eval_type
    {
    };

    /**
     * @brief Tag to indicate that derivative of the spline should be evaluated.
     */
    struct eval_deriv_type
    {
    };

public:
    /**
     * @brief Tag the B-spline type of the first dimension.
     */
    using bsplines_type1 = BSplinesType1;
    /**
     * @brief Tag the B-spline type of the second dimension.
     */
    using bsplines_type2 = BSplinesType2;
    /**
     * @brief Tag the first dimension.
     */
    using Dim1 = typename BSplinesType1::tag_type;
    /**
     * @brief Tag the second dimension.
     */
    using Dim2 = typename BSplinesType2::tag_type;

private:
    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_left_bc_1;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_right_bc_1;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_left_bc_2;

    SplineBoundaryValue2D<BSplinesType1, BSplinesType2> const& m_right_bc_2;

public:
    SplineEvaluator2D() = delete;

    /**
     * @brief Instantiate an evaluator operator.
     *
     * @param[in] left_bc_1
     * 			A SplineBoundaryValue2D object giving the value on the "left side" of the domain
     * 			in the first dimension.
     * @param[in] right_bc_1
     * 			A SplineBoundaryValue2D object giving the value on the "right side" of the domain
     * 			in the first dimension.
     * @param[in] left_bc_2
     * 			A SplineBoundaryValue2D object giving the value on the "left side" of the domain
     * 			in the second dimension.
     * @param[in] right_bc_2
     * 			A SplineBoundaryValue2D object giving the value on the "right side" of the domain
     * 			in the second dimension.
     *
     * @see SplineBoundaryValue2D
     */
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

    /**
     * @brief Instantiate a SplineEvaluator2D from another
     * SplineEvaluator2D (lvalue).
     *
     * @param[in] x
     * 		SplineEvaluator2D evaluator used to instantiate the new one.
     */
    SplineEvaluator2D(SplineEvaluator2D const& x) = default;

    /**
     * @brief Instantiate a SplineEvaluator2D from another temporary
     * SplineEvaluator2D (rvalue).
     *
     * @param[in] x
     * 		SplineEvaluator2D evaluator used to instantiate the new one.
     */
    SplineEvaluator2D(SplineEvaluator2D&& x) = default;

    ~SplineEvaluator2D() = default;

    /**
     * @brief Assign a SplineEvaluator2D from another SplineEvaluator2D (lvalue).
     *
     * @param[in] x
     * 		SplineEvaluator2D mapping used to assign.
     *
     * @return The SplineEvaluator2D assigned.
     */
    SplineEvaluator2D& operator=(SplineEvaluator2D const& x) = default;

    /**
     * @brief Assign a SplineEvaluator2D from another temporary SplineEvaluator2D (rvalue).
     *
     * @param[in] x
     * 		SplineEvaluator2D mapping used to assign.
     *
     * @return The SplineEvaluator2D assigned.
     */
    SplineEvaluator2D& operator=(SplineEvaluator2D&& x) = default;

    /**
     * @brief Get the value of the function on B-splines at the coordinate given.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     *
     * @return A double containing the value of the function at the coordinate given.
     */
    double operator()(
            ddc::Coordinate<Dim1, Dim2> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        return eval(coord_eval, spline_coef, vals1, vals2);
    }

    /**
     * @brief Get the values of the function on B-splines at the coordinates given.
     *
     * @param[out] spline_eval
     * 			A ChunkSpan containing the values of the function evaluated at the coordinates given.
     * @param[in] coords_eval
     * 			A ChunkSpan with the 2D coordinates where we want to evaluate the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     */
    template <class Domain>
    void operator()(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval(coords_eval(i), spline_coef, vals1, vals2);
        });
    }

    /**
     * @brief Get the value of the derivative of the first dimension of the function on B-splines at the coordinate given.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate the derivative of the first dimension of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     *
     * @return A double containing the value of the derivative of the first dimension of the function at the coordinate given.
     */
    double deriv_dim_1(
            ddc::Coordinate<Dim1, Dim2> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        return eval_no_bc(
                ddc::select<Dim1>(coord_eval),
                ddc::select<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_deriv_type(),
                eval_type());
    }

    /**
     * @brief Get the value of the derivative of the second dimension of the function on B-splines at the coordinate given.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate the derivative of the second dimension of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     *
     * @return A double containing the value of the derivative of the second dimension of the function at the coordinate given.
     */
    double deriv_dim_2(
            ddc::Coordinate<Dim1, Dim2> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        return eval_no_bc(
                ddc::select<Dim1>(coord_eval),
                ddc::select<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_type(),
                eval_deriv_type());
    }

    /**
     * @brief Get the value of the cross derivative of the function on B-splines at the coordinate given.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate the cross derivative of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     *
     * @return A double containing the value of the cross derivative of the function at the coordinate given.
     */
    double deriv_1_and_2(
            ddc::Coordinate<Dim1, Dim2> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        return eval_no_bc(
                ddc::select<Dim1>(coord_eval),
                ddc::select<Dim2>(coord_eval),
                spline_coef,
                vals1,
                vals2,
                eval_deriv_type(),
                eval_deriv_type());
    }

    /**
     * @brief Get the values of the derivative of the first dimension of the function on B-splines at the coordinates given.
     *
     * @param[out] spline_eval
     * 			A ChunkSpan with the values of the derivative of the first dimension of the function at the coordinates given.
     * @param[in] coords_eval
     * 			A ChunkSpan with the 2D coordinates where we want to evaluate the derivative of the first dimension of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     */
    template <class Domain>
    void deriv_dim_1(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    ddc::select<Dim1>(coords_eval(i)),
                    ddc::select<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_deriv_type(),
                    eval_type());
        });
    }

    /**
     * @brief Get the values of the derivative of the second dimension of the function on B-splines at the coordinates given.
     *
     * @param[out] spline_eval
     * 			A ChunkSpan with the values of the derivative of the second dimension of the function at the coordinates given.
     * @param[in] coords_eval
     * 			A ChunkSpan with the 2D coordinates where we want to evaluate the derivative of the second dimension of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     */
    template <class Domain>
    void deriv_dim_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    ddc::select<Dim1>(coords_eval(i)),
                    ddc::select<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_type(),
                    eval_deriv_type());
        });
    }

    /**
     * @brief Get the values of the cross derivative of the function on B-splines at the coordinates given.
     *
     * @param[out] spline_eval
     * 			A ChunkSpan with the values of the cross derivative of the function at the coordinates given.
     * @param[in] coords_eval
     * 			A ChunkSpan with the 2D coordinates where we want to evaluate the cross derivative of the function.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     */
    template <class Domain>
    void deriv_dim_1_and_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<Dim1, Dim2> const, Domain> const coords_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        std::array<double, bsplines_type1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, bsplines_type2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(
                    ddc::select<Dim1>(coords_eval(i)),
                    ddc::select<Dim2>(coords_eval(i)),
                    spline_coef,
                    vals1,
                    vals2,
                    eval_deriv_type(),
                    eval_deriv_type());
        });
    }

    /**
     * @brief Get the the integral of the function on B-splines on the domain.
     *
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to integrate.
     *
     * @return A double with the value of the integral of the function.
     */
    double integrate(
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef) const
    {
        ddc::Chunk<double, ddc::DiscreteDomain<BSplinesType1>> values1(
                ddc::DiscreteDomain<BSplinesType1>(spline_coef.domain()));
        DSpan1D vals1 = values1.allocation_mdspan();
        ddc::Chunk<double, ddc::DiscreteDomain<BSplinesType2>> values2(
                ddc::DiscreteDomain<BSplinesType2>(spline_coef.domain()));
        DSpan1D vals2 = values2.allocation_mdspan();

        ddc::discrete_space<bsplines_type1>().integrals(values1.span_view());
        ddc::discrete_space<bsplines_type2>().integrals(values2.span_view());

        return ddc::transform_reduce(
                spline_coef.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                [&](ddc::DiscreteElement<BSplinesType1, BSplinesType2> const i) {
                    return spline_coef(i) * values1(ddc::select<BSplinesType1>(i))
                           * values2(ddc::select<BSplinesType2>(i));
                });
    }

private:
    /**
     * @brief Evaluate the function on B-splines at the coordinate given.
     *
     * This function firstly deals with the boundary conditions and calls the SplineEvaluator2D::eval_no_bc function
     * to evaluate.
     *
     * @param[in] coord_eval
     * 			The 2D coordinate where we want to evaluate.
     * @param[in] spline_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     * @param[out] vals1
     * 			A ChunkSpan with the not-null values of each function of the spline in the first dimension.
     * @param[out] vals2
     * 			A ChunkSpan with the not-null values of each function of the spline in the second dimension.
     *
     * @return A double with the value of the function at the coordinate given.
     *
     * @see SplineBoundaryValue
     */
    double eval(
            ddc::Coordinate<Dim1, Dim2> const& coord_eval,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef,
            DSpan1D const vals1,
            DSpan1D const vals2) const
    {
        ddc::Coordinate<Dim1> coord_eval1 = ddc::select<Dim1>(coord_eval);
        ddc::Coordinate<Dim2> coord_eval2 = ddc::select<Dim2>(coord_eval);
        if constexpr (bsplines_type1::is_periodic()) {
            if (coord_eval1 < ddc::discrete_space<bsplines_type1>().rmin()
                || coord_eval1 > ddc::discrete_space<bsplines_type1>().rmax()) {
                coord_eval1 -= std::floor(
                                       (coord_eval1 - ddc::discrete_space<bsplines_type1>().rmin())
                                       / ddc::discrete_space<bsplines_type1>().length())
                               * ddc::discrete_space<bsplines_type1>().length();
            }
        }

        if constexpr (bsplines_type2::is_periodic()) {
            if (coord_eval2 < ddc::discrete_space<bsplines_type2>().rmin()
                || coord_eval2 > ddc::discrete_space<bsplines_type2>().rmax()) {
                coord_eval2 -= std::floor(
                                       (coord_eval2 - ddc::discrete_space<bsplines_type2>().rmin())
                                       / ddc::discrete_space<bsplines_type2>().length())
                               * ddc::discrete_space<bsplines_type2>().length();
            }
        }

        if constexpr (!bsplines_type1::is_periodic()) {
            if (coord_eval1 < ddc::discrete_space<bsplines_type1>().rmin()) {
                return m_left_bc_1(coord_eval1, coord_eval2, spline_coef);
            }
            if (coord_eval1 > ddc::discrete_space<bsplines_type1>().rmax()) {
                return m_right_bc_1(coord_eval1, coord_eval2, spline_coef);
            }
        }

        if constexpr (!bsplines_type2::is_periodic()) {
            if (coord_eval2 < ddc::discrete_space<bsplines_type2>().rmin()) {
                return m_left_bc_2(coord_eval1, coord_eval2, spline_coef);
            }
            if (coord_eval2 > ddc::discrete_space<bsplines_type2>().rmax()) {
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

    /**
     * @brief Evaluate the function or its derivative at the coordinate given.
     *
     * @param[in] coord_eval1
     * 			The coordinate on the first dimension where we want to evaluate.
     * @param[in] coord_eval2
     * 			The coordinate on the second dimension where we want to evaluate.
     * @param[in] splne_coef
     * 			The B-splines coefficients of the function we want to evaluate.
     * @param[out] vals1
     * 			A ChunkSpan with the not-null values of each function of the spline in the first dimension.
     * @param[out] vals2
     * 			A ChunkSpan with the not-null values of each function of the spline in the second dimension.
     * @param[in] eval_type_1
     * 			A flag indicating if we evaluate the function or its derivative in the first dimension.
     * 			The type of this object is either `eval_type` or `eval_deriv_type`.
     * @param[in] eval_type_2
     * 			A flag indicating if we evaluate the function or its derivative in the second dimension.
     *          The type of this object is either `eval_type` or `eval_deriv_type`.
     */
    template <class EvalType1, class EvalType2>
    double eval_no_bc(
            ddc::Coordinate<Dim1> const& coord_eval1,
            ddc::Coordinate<Dim2> const& coord_eval2,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplinesType1, BSplinesType2>> const
                    spline_coef,
            DSpan1D const vals1,
            DSpan1D const vals2,
            EvalType1 const eval_type_1,
            EvalType2 const eval_type_2) const
    {
        static_assert(
                std::is_same_v<EvalType1, eval_type> || std::is_same_v<EvalType1, eval_deriv_type>);
        static_assert(
                std::is_same_v<EvalType2, eval_type> || std::is_same_v<EvalType2, eval_deriv_type>);
        ddc::DiscreteElement<BSplinesType1> jmin1;
        ddc::DiscreteElement<BSplinesType2> jmin2;

        if constexpr (std::is_same_v<EvalType1, eval_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>().eval_basis(vals1, coord_eval1);
        } else if constexpr (std::is_same_v<EvalType1, eval_deriv_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>().eval_deriv(vals1, coord_eval1);
        }
        if constexpr (std::is_same_v<EvalType2, eval_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>().eval_basis(vals2, coord_eval2);
        } else if constexpr (std::is_same_v<EvalType2, eval_deriv_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>().eval_deriv(vals2, coord_eval2);
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
