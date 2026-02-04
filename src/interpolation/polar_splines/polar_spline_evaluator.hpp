// SPDX-License-Identifier: MIT
#pragma once

#include "view.hpp"

/**
 * @brief Define an evaluator on polar B-splines.
 *
 * @see PolarBSplines
 */
template <class ExecSpace, class MemorySpace, class PolarBSplinesType, class OuterExtrapolationRule>
class PolarSplineEvaluator
{
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
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;
    /**
     * @brief Tag the first dimension of the space.
     */
    using R = typename BSplinesR::continuous_dimension_type;
    /**
     * @brief Tag the second dimension of the space.
     */
    using Theta = typename BSplinesTheta::continuous_dimension_type;

    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

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
    KOKKOS_FUNCTION double operator()(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
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
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval(coords_eval(i), spline_coef);
                });
    }

    /**
     * @brief Get the values of the spline function on a domain.
     *
     * @param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate.
     */
    template <class Domain>
    void operator()(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval(ddc::coordinate(i), spline_coef);
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
     * @param[in] deriv_order
     *      The index of the derivative order (e.g. Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1,3)
     *      for the cross derivative @f$ dr d \theta^3 @f$.
     *
     * @return The value of the derivative of the spline function on the
     * first dimension.
     */
    template <class DerivDim>
    KOKKOS_FUNCTION double deriv(
            Idx<ddc::Deriv<DerivDim>> deriv_order,
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, deriv_order);
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
    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] KOKKOS_FUNCTION double
    deriv_dim_1(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<R>>(1));
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
    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] KOKKOS_FUNCTION double
    deriv_dim_2(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<Theta>>(1));
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
    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) "
                 "instead")]] KOKKOS_FUNCTION double
    deriv_1_and_2(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
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
     * @param[in] deriv_order
     *      The index of the derivative order (e.g. Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1,3)
     *      for the cross derivative @f$ dr d \theta^3 @f$.
     */
    template <class Domain, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const deriv_order,
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, deriv_order);
                });
    }

    /**
     * @brief Get the values of the derivative of the spline function on the
     * first dimension.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate.
     * @param[in] deriv_order
     *      The index of the derivative order (e.g. Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1,3)
     *      for the cross derivative @f$ dr d \theta^3 @f$.
     */
    template <class Domain, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const deriv_order,
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(ddc::coordinate(i), spline_coef, deriv_order);
                });
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
    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] void deriv_dim_1(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, Idx<ddc::Deriv<R>>(1));
                });
    }

    /**
     * @brief Get the values of the derivative of the spline function on the
     * first dimension.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate.
     */
    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] void deriv_dim_1(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i)
                            = eval_no_bc(ddc::coordinate(i), spline_coef, Idx<ddc::Deriv<R>>(1));
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
    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] void deriv_dim_2(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i)
                            = eval_no_bc(coords_eval(i), spline_coef, Idx<ddc::Deriv<Theta>>(1));
                });
    }

    /**
     * @brief Get the values of the derivative of the spline function on the
     * second dimension.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] spline_coef
     *      The B-splines coefficients of the spline function we want to evaluate..
     */
    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] void deriv_dim_2(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            ddc::coordinate(i),
                            spline_coef,
                            Idx<ddc::Deriv<Theta>>(1));
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
    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) instead")]] void
    deriv_dim_1_and_2(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            coords_eval(i),
                            spline_coef,
                            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
                });
    }

    /**
     * @brief Get the values of the cross derivative of the spline function.
     *
     *@param[out] spline_eval
     *      The values of the function evaluated on the domain.
     * @param[in] spline_coef
     *      The B-splines coefficients of the splinefunction we want to evaluate.
     */
    template <class Domain>
    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) instead")]] void
    deriv_dim_1_and_2(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            ddc::coordinate(i),
                            spline_coef,
                            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
                });
    }

private:
    KOKKOS_FUNCTION double eval(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        const double coord_eval1 = ddc::get<R>(coord_eval);
        double coord_eval2 = ddc::get<Theta>(coord_eval);
        if (coord_eval1 > ddc::discrete_space<BSplinesR>().rmax()) {
            return m_outer_bc(coord_eval, spline_coef);
        }
        if (coord_eval2 < ddc::discrete_space<BSplinesTheta>().rmin()
            || coord_eval2 > ddc::discrete_space<BSplinesTheta>().rmax()) {
            coord_eval2 -= Kokkos::floor(
                                   (coord_eval2 - ddc::discrete_space<BSplinesTheta>().rmin())
                                   / ddc::discrete_space<BSplinesTheta>().length())
                           * ddc::discrete_space<BSplinesTheta>().length();
        }
        Coord<R, Theta> coord_eval_new(coord_eval1, coord_eval2);
        return eval_no_bc(coord_eval_new, spline_coef, Idx<>());
    }

    template <class... DerivDims>
    KOKKOS_FUNCTION double eval_no_bc(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef,
            Idx<DerivDims...> const deriv_order) const
    {
        std::array<double, PolarBSplinesType::n_singular_basis()> singular_data;
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, (BSplinesR::degree() + 1) * (BSplinesTheta::degree() + 1)> data;
        DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesTheta::degree() + 1);

        Idx<BSplinesR, BSplinesTheta> jmin
                = ddc::discrete_space<PolarBSplinesType>()
                          .eval_deriv(singular_vals, vals, coord_eval, deriv_order);

        double y = 0.0;
        for (std::size_t i = 0; i < PolarBSplinesType::n_singular_basis(); ++i) {
            y += spline_coef(Idx<PolarBSplinesType>(i)) * singular_vals(i);
        }
        Idx<BSplinesR> jmin_r(jmin);
        Idx<BSplinesTheta> jmin_theta(jmin);

        DConstField<IdxRange<BSplinesR, BSplinesTheta>, MemorySpace> spline_coef_2d
                = PolarBSplinesType::get_tensor_product_subset(spline_coef);
        IdxRange<BSplinesR, BSplinesTheta> tensor_prod_idx_range = get_idx_range(spline_coef_2d);
        IdxRange<BSplinesTheta> tensor_prod_idx_range_theta(tensor_prod_idx_range);
        Idx<BSplinesTheta> idx_theta_max = tensor_prod_idx_range_theta.back();
        IdxStep<BSplinesTheta> n_idx_theta = tensor_prod_idx_range_theta.extents();
        for (int i = 0; i < BSplinesR::degree() + 1; ++i) {
            for (std::size_t j = 0; j < BSplinesTheta::degree() + 1; ++j) {
                Idx<BSplinesTheta> idx_theta = jmin_theta + j;
                if (idx_theta > idx_theta_max) {
                    idx_theta -= n_idx_theta;
                }
                y += spline_coef_2d(jmin_r + i, idx_theta) * vals(i, j);
            }
        }
        return y;
    }
};
