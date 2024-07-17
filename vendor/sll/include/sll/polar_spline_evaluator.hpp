#pragma once

#include <sll/polar_spline.hpp>
#include <sll/view.hpp>

/**
 * @brief Define an evaluator on polar B-splines.
 *
 * @see PolarBSplines
 */
template <class PolarBSplinesType, class OuterExtrapolationRule>
class PolarSplineEvaluator
{
private:
    // Tags to determine what to evaluate
    /**
     * @brief Tag for the evaluation of the function.
     */
    struct eval_type
    {
    };

    /**
     * @brief Tag for the evaluation of the derivative on the
     * first dimension.
     */
    struct eval_deriv_r_type
    {
    };

    /**
     * @brief Tag for the evaluation of the derivative on the
     * second dimension.
     */
    struct eval_deriv_p_type
    {
    };

    /**
     * @brief Tag for the evaluation of the cross derivative of the function.
     */
    struct eval_deriv_r_p_type
    {
    };

public:
    /**
     * @brief Tag the type of B-splines.
     */
    using bsplines_type = PolarBSplinesType;
    /**
     * @brief Tag the type of first dimension B-splines.
     */
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    /**
     * @brief Tag the type of second dimension B-splines.
     */
    using BSplinesP = typename PolarBSplinesType::BSplinesP_tag;
    /**
     * @brief Tag the first dimension of the space.
     */
    using DimR = typename BSplinesR::continuous_dimension_type;
    /**
     * @brief Tag the second dimension of the space.
     */
    using DimP = typename BSplinesP::continuous_dimension_type;

public:
    /**
     * @brief Tag the order of continuity of the B-splines basis
     * at the O-point.
     */
    static int constexpr continuity = PolarBSplinesType::continuity;

private:
    OuterExtrapolationRule m_outer_bc;

public:
    PolarSplineEvaluator() = delete;

    /**
     * @brief Instantiate a PolarSplineEvaluator with boundary evaluation conditions.
     *
     * Instantiate a PolarSplineEvaluator by specifying how points lying outside the
     * domain should be evaluated. The domain is the domain on which the polar splines
     * are defined.
     *
     * @param[in] outer_bc
     *      A class containing an operator which can be called to provide a boundary value to
     *      evaluate a point lying outside the domain.
     */
    explicit PolarSplineEvaluator(OuterExtrapolationRule const& outer_bc) : m_outer_bc(outer_bc) {}

    /**
     * @brief Instantiate a PolarSplineEvaluator from another.
     *
     * @param[in] x
     *      Another PolarSplineEvaluator (lvalue).
     */
    PolarSplineEvaluator(PolarSplineEvaluator const& x) = default;

    /**
     * @brief Instantiate a PolarSplineEvaluator from another temporary.
     *
     * @param[in] x
     *      Another temporary PolarSplineEvaluator (rvalue).
     */
    PolarSplineEvaluator(PolarSplineEvaluator&& x) = default;

    ~PolarSplineEvaluator() = default;

    /**
     * @brief Assign a PolarSplineEvaluator from another.
     *
     * @param[in] x
     *      Another PolarSplineEvaluator (lvalue).
     *
     * @return PolarSplineEvaluator assigned.
     */
    PolarSplineEvaluator& operator=(PolarSplineEvaluator const& x) = default;

    /**
     * @brief Assign a PolarSplineEvaluator from another temporary.
     *
     * @param[in] x
     *      Another temporary PolarSplineEvaluator (rvalue).
     *
     * @return PolarSplineEvaluator assigned.
     */
    PolarSplineEvaluator& operator=(PolarSplineEvaluator&& x) = default;

    /**
     * @brief Get the value of the spline function at a given coordinate.
     *
     * @param[in] coord_eval
     *      The coordinate where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the function we want to evaluate.
     *
     * @return A double with value of the spline function at the given coordinate.
     */
    double operator()(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    /**
     * @brief Get the values of the spline function on a domain.
     *
     * @param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] coords_eval
     *      The coordinates where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate.
     */
    template <class Domain>
    void operator()(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval(coords_eval(i), spline_coef);
        });
    }

    /**
     * @brief Get the value of the derivative of the spline function on the
     * first dimension.
     *
     * @param[in] coord_eval
     *      The coordinate where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the function we want to evaluate.
     *
     * @return The value of the derivative of the spline function on the
     * first dimension.
     */
    double deriv_dim_1(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_type());
    }

    /**
     * @brief Get the value of the derivative of the spline function on the
     * second dimension.
     *
     * @param[in] coord_eval
     *      The coordinate where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the function we want to evaluate.
     *
     * @return The value of the derivative of the spline function on the
     * second dimension.
     */
    double deriv_dim_2(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_p_type());
    }


    /**
     * @brief Get the value of the cross derivative of the spline function.
     *
     * @param[in] coord_eval
     *      The coordinate where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the function we want to evaluate.
     *
     * @return The value of the cross derivative of the spline
     * function
     */
    double deriv_1_and_2(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_p_type());
    }

    /**
     * @brief Get the values of the derivative of the spline function on the
     * first dimension.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] coords_eval
     *      The coordinates where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate.
     */
    template <class Domain>
    void deriv_dim_1(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_type());
        });
    }

    /**
     * @brief Get the values of the derivative of the spline function on the
     * second dimension.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] coords_eval
     *      The coordinates where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate..
     */
    template <class Domain>
    void deriv_dim_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_p_type());
        });
    }

    /**
     * @brief Get the values of the cross derivative of the spline function.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] coords_eval
     *      The coordinates where we want to evaluate.
     * @param[in] spline_coef
     *      The B-splines coefficients of the splinefunction we want to evaluate.
     */
    template <class Domain>
    void deriv_dim_1_and_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_p_type());
        });
    }

    /**
     * @brief Get the in integral of a spline function over the domain.
     *
     * @param[in] spline_coef
     *      The B-splines coefficients of the function we want to evaluate.
     * @param[in] mapping
     *      The mapping function.
     *
     * @return The integral of the spline function over the domain.
     */
    template <class Mapping>
    double integrate(PolarSplineView<PolarBSplinesType> const spline_coef, Mapping const mapping)
            const
    {
        int constexpr nr = ddc::discrete_space<BSplinesR>().ncells() + BSplinesR::degree() - 2;
        int constexpr np = ddc::discrete_space<BSplinesP>().ncells() + BSplinesP::degree();
        std::array<double, PolarBSplinesType::eval_size()> singular_values;
        DSpan1D singular_vals(singular_values.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, nr * np> values;
        DSpan2D vals(values.data(), nr, np);

        ddc::discrete_space<PolarBSplinesType>().integrals(singular_vals, vals);

        double y = 0.;
        ddc::for_each(
                spline_coef.singular_spline_coef.domain(),
                [=](ddc::DiscreteElement<PolarBSplinesType> const i) {
                    y += spline_coef.singular_spline_coef(i) * singular_vals(i)
                         * mapping.determinant(i);
                });
        ddc::for_each(
                spline_coef.spline_coef.domain(),
                [=](ddc::DiscreteElement<BSplinesR, BSplinesP> const i) {
                    y += spline_coef.spline_coef(i)
                         * vals(ddc::select<BSplinesR>(i), ddc::select<BSplinesP>(i))
                         * mapping.determinant(i);
                });
        return y;
    }

private:
    double eval(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        const double coord_eval1 = ddc::get<DimR>(coord_eval);
        double coord_eval2 = ddc::get<DimP>(coord_eval);
        if (coord_eval1 > ddc::discrete_space<BSplinesR>().rmax()) {
            return m_outer_bc(coord_eval, spline_coef);
        }
        if (coord_eval2 < ddc::discrete_space<BSplinesP>().rmin()
            || coord_eval2 > ddc::discrete_space<BSplinesP>().rmax()) {
            coord_eval2 -= std::floor(
                                   (coord_eval2 - ddc::discrete_space<BSplinesP>().rmin())
                                   / ddc::discrete_space<BSplinesP>().length())
                           * ddc::discrete_space<BSplinesP>().length();
        }
        return eval_no_bc(coord_eval, spline_coef, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef,
            EvalType const) const
    {
        static_assert(
                std::is_same_v<
                        EvalType,
                        eval_type> || std::is_same_v<EvalType, eval_deriv_r_type> || std::is_same_v<EvalType, eval_deriv_p_type> || std::is_same_v<EvalType, eval_deriv_r_p_type>);

        std::array<double, PolarBSplinesType::n_singular_basis()> singular_data;
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, (BSplinesR::degree() + 1) * (BSplinesP::degree() + 1)> data;
        DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesP::degree() + 1);

        ddc::DiscreteElement<BSplinesR, BSplinesP> jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_basis(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_p_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_p(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_p_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r_and_p(singular_vals, vals, coord_eval);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < PolarBSplinesType::n_singular_basis(); ++i) {
            y += spline_coef.singular_spline_coef(ddc::DiscreteElement<PolarBSplinesType>(i))
                 * singular_vals(i);
        }
        ddc::DiscreteElement<BSplinesR> jmin_r = ddc::select<BSplinesR>(jmin);
        ddc::DiscreteElement<BSplinesP> jmin_p = ddc::select<BSplinesP>(jmin);
        int nr = BSplinesR::degree() + 1;
        if (jmin_r.uid() < continuity + 1) {
            nr = nr - (continuity + 1 - jmin_r.uid());
            jmin_r = ddc::DiscreteElement<BSplinesR>(continuity + 1);
        }
        for (int i = 0; i < nr; ++i) {
            for (std::size_t j = 0; j < BSplinesP::degree() + 1; ++j) {
                y += spline_coef.spline_coef(jmin_r + i, jmin_p + j) * vals(i, j);
            }
        }
        return y;
    }
};
