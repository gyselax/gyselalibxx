// SPDX-License-Identifier: MIT
#pragma once

/**
 * @brief A class to evaluate, differentiate or integrate a Lagrange function.
 *
 * A class which contains an operator () which can be used to evaluate, or
 * differentiate a Lagrange polynomial.
 *
 * @tparam ExecSpace The Kokkos execution space on which the Lagrange evaluation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (Lagrange coefficients and evaluation) is stored.
 * @tparam DataType The data type on which calculations are made.
 * @tparam LagrangeBasis The discrete dimension representing the Lagrange basis.
 * @tparam InterpolationGrid The discrete dimension on which evaluation points are defined.
 * @tparam LowerExtrapolationRule The lower extrapolation rule type.
 * @tparam UpperExtrapolationRule The upper extrapolation rule type.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class LagrangeBasis,
        class InterpolationGrid,
        class LowerExtrapolationRule,
        class UpperExtrapolationRule>
class LagrangeEvaluator
{
public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The type of the evaluation continuous dimension (continuous dimension of interest) used by this class.
    using continuous_dimension_type = typename LagrangeBasis::continuous_dimension_type;

    /// @brief The discrete dimension representing the Lagrange basis.
    using lagrange_basis_type = LagrangeBasis;

    using basis_domain_type =
            typename LagrangeBasis::template Impl<LagrangeBasis, MemorySpace>::knot_grid;

    /// @brief The type of the domain for the 1D evaluation mesh used by this class.
    using evaluation_domain_type = IdxRange<InterpolationGrid>;

    /**
     * @brief The type of the whole domain representing evaluation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_evaluation_domain_type = BatchedInterpolationGrid;

    /// @brief The type of the 1D Lagrange domain corresponding to the dimension of interest.
    using lagrange_domain_type = IdxRange<basis_domain_type>;

    /**
     * @brief The type of the batch domain (obtained by removing the dimension of interest
     * from the whole domain).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batch_domain_type = ddc::remove_dims_of_t<BatchedInterpolationGrid, InterpolationGrid>;

    /**
     * @brief The type of the whole Lagrange domain (cartesian product of 1D Lagrange domain
     * and batch domain) preserving the order of dimensions.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationDDom,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationDDom>>>
    using batched_lagrange_domain_type
            = ddc::replace_dim_of_t<BatchedInterpolationDDom, InterpolationGrid, basis_domain_type>;

    /// @brief The type of the extrapolation rule at the lower boundary.
    using lower_extrapolation_rule_type = LowerExtrapolationRule;

    /// @brief The type of the extrapolation rule at the upper boundary.
    using upper_extrapolation_rule_type = UpperExtrapolationRule;


private:
    LowerExtrapolationRule m_lower_extrap_rule;

    UpperExtrapolationRule m_upper_extrap_rule;

public:
    static_assert(
            std::is_same_v<LowerExtrapolationRule,
                            typename ddc::PeriodicExtrapolationRule<
                                    continuous_dimension_type>> == lagrange_basis_type::is_periodic()
                    && std::is_same_v<
                               UpperExtrapolationRule,
                               typename ddc::PeriodicExtrapolationRule<
                                       continuous_dimension_type>> == lagrange_basis_type::is_periodic(),
            "PeriodicExtrapolationRule has to be used if and only if dimension is periodic");
    static_assert(
            std::is_invocable_r_v<
                    DataType,
                    LowerExtrapolationRule,
                    Coord<continuous_dimension_type>,
                    ConstField<DataType, lagrange_domain_type, memory_space>>,
            "LowerExtrapolationRule::operator() has to be callable with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    DataType,
                    UpperExtrapolationRule,
                    Coord<continuous_dimension_type>,
                    ConstField<DataType, lagrange_domain_type, memory_space>>,
            "UpperExtrapolationRule::operator() has to be callable with usual arguments.");

    /**
     * @brief Build a LagrangeEvaluator acting on batched_lagrange_domain.
     *
     * @param lower_extrap_rule The extrapolation rule at the lower boundary.
     * @param upper_extrap_rule The extrapolation rule at the upper boundary.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    explicit LagrangeEvaluator(
            LowerExtrapolationRule const& lower_extrap_rule,
            UpperExtrapolationRule const& upper_extrap_rule)
        : m_lower_extrap_rule(lower_extrap_rule)
        , m_upper_extrap_rule(upper_extrap_rule)
    {
    }

    /**
     * @brief Copy-constructs.
     *
     * @param x A reference to another LagrangeEvaluator.
     */
    LagrangeEvaluator(LagrangeEvaluator const& x) = default;

    /**
     * @brief Move-constructs.
     *
     * @param x An rvalue to another LagrangeEvaluator.
     */
    LagrangeEvaluator(LagrangeEvaluator&& x) = default;

    /// @brief Destructs
    ~LagrangeEvaluator() = default;

    /**
     * @brief Copy-assigns.
     *
     * @param x A reference to another LagrangeEvaluator.
     * @return A reference to this object.
     */
    LagrangeEvaluator& operator=(LagrangeEvaluator const& x) = default;

    /**
     * @brief Move-assigns.
     *
     * @param x An rvalue to another LagrangeEvaluator.
     * @return A reference to this object.
     */
    LagrangeEvaluator& operator=(LagrangeEvaluator&& x) = default;

    /**
     * @brief Get the lower extrapolation rule.
     *
     * Extrapolation rules are functors used to define the behaviour of the LagrangeEvaluator
     * outside the domain where the break points of the Lagrange basis are defined.
     *
     * @return The lower extrapolation rule.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    lower_extrapolation_rule_type lower_extrapolation_rule() const
    {
        return m_lower_extrap_rule;
    }

    /**
     * @brief Get the upper extrapolation rule.
     *
     * Extrapolation rules are functors used to define the behaviour of the LagrangeEvaluator
     * outside the domain where the break points of the Lagrange basis are defined.
     *
     * @return The upper extrapolation rule.
     *
     * @see NullExtrapolationRule ConstantExtrapolationRule PeriodicExtrapolationRule
     */
    upper_extrapolation_rule_type upper_extrapolation_rule() const
    {
        return m_upper_extrap_rule;
    }

    /**
     * @brief Evaluate 1D Lagrange polynomial (described by its Lagrange coefficients) at a given coordinate.
     *
     * The Lagrange coefficients describe 1D Lagrange polynomials defined on a Lagrange basis.
     * They can be obtained using the IdentityInterpolationBuilder.
     *
     * @param coord_eval The coordinate where the Lagrange polynomial is evaluated.
     *              Note that only the component along the dimension of interest is used.
     * @param lagrange_coef A Field storing the 1D Lagrange coefficients.
     *
     * @return The value of the Lagrange polynomial at the desired coordinate.
     */
    template <class Layout, class... CoordsDims, class BatchedLagrangeIdxRange>
    KOKKOS_FUNCTION DataType operator()(
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, BatchedLagrangeIdxRange, memory_space, Layout> const lagrange_coef)
            const
    {
        return eval(coord_eval, lagrange_coef);
    }

    /**
     * @brief Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh.
     *
     * The Lagrange coefficients describe Lagrange polynomials defined on a cartesian product of batch_idx_range and
     * Lagrange basis. They can be obtained using the IdentityInterpolationBuilder.
     *
     * This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates
     * identified by a batch_domain_type::discrete_element_type, the evaluation is performed with the 1D set of
     * Lagrange coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @param[out] lagrange_eval The values of the Lagrange polynomials at the desired coordinates.
     * @param[in] coords_eval The coordinates where the Lagrange polynomials are evaluated. Those are
     * stored in a Field defined on a BatchedInterpolationIdxRange. Note that the coordinates of the
     * points represented by this index range are unused and irrelevant.
     * @param[in] lagrange_coef A Field storing the Lagrange coefficients.
     */
    template <
            class Layout1,
            class Layout2,
            class Layout3,
            class BatchedInterpolationIdxRange,
            class... CoordsDims>
    void operator()(
            Field<DataType, BatchedInterpolationIdxRange, memory_space, Layout1> const
                    lagrange_eval,
            ConstField<
                    Coord<CoordsDims...>,
                    BatchedInterpolationIdxRange,
                    memory_space,
                    Layout2> const coords_eval,
            ConstField<
                    DataType,
                    batched_lagrange_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout3> const lagrange_coef) const
    {
        evaluation_domain_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_domain_type<BatchedInterpolationGrid> const batch_idx_range(
                get_idx_range(lagrange_eval));

        using BatchedInterpolationIdx =
                typename batch_domain_type<BatchedInterpolationGrid>::discrete_element_type;

        ddc::parallel_for_each(
                "Lagrange_evaluate",
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(BatchedInterpolationIdx const j) {
                    // auto allows non-contiguous layouts
                    auto const lagrange_eval_1D = lagrange_eval[j];
                    auto const coords_eval_1D = coords_eval[j];
                    auto const lagrange_coef_1D = lagrange_coef[j];
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        lagrange_eval_1D(i) = eval(coords_eval_1D(i), lagrange_coef_1D);
                    }
                });
    }

    /**
     * @brief Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh.
     *
     * The Lagrange coefficients describe Lagrange polynomials defined on a cartesian product of batch_idx_range and
     * Lagrange basis. They can be obtained using the IdentityInterpolationBuilder.
     *
     * This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates
     * identified by a batch_domain_type::discrete_element_type, the evaluation is performed with the 1D set of
     * Lagrange coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * @param[out] lagrange_eval The values of the Lagrange polynomials at the mesh points.
     * @param[in] lagrange_coef A Field storing the Lagrange coefficients.
     */
    template <class Layout1, class Layout2, class BatchedInterpolationGrid>
    void operator()(
            Field<DataType, BatchedInterpolationGrid, memory_space, Layout1> const lagrange_eval,
            ConstField<
                    DataType,
                    batched_lagrange_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout2> const lagrange_coef) const
    {
        evaluation_domain_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_domain_type<BatchedInterpolationGrid> const batch_idx_range(
                get_idx_range(lagrange_eval));

        ddc::parallel_for_each(
                "lagrange_evaluate",
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type<
                                    BatchedInterpolationGrid>::discrete_element_type const j) {
                    auto const lagrange_eval_1D = lagrange_eval[j];
                    auto const lagrange_coef_1D = lagrange_coef[j];
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        Coord<continuous_dimension_type> coord_eval_1D = ddc::coordinate(i);
                        lagrange_eval_1D(i) = eval(coord_eval_1D, lagrange_coef_1D);
                    }
                });
    }

private:
    template <class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType
    eval(Coord<CoordsDims...> const& coord_eval,
         ConstField<DataType, lagrange_domain_type, memory_space, Layout> const lagrange_coef) const
    {
        Coord<continuous_dimension_type> coord_eval_interest(coord_eval);
        if constexpr (lagrange_basis_type::is_periodic()) {
            if (coord_eval_interest < ddc::discrete_space<lagrange_basis_type>().rmin()
                || coord_eval_interest > ddc::discrete_space<lagrange_basis_type>().rmax()) {
                coord_eval_interest
                        -= Kokkos::floor(
                                   (coord_eval_interest
                                    - ddc::discrete_space<lagrange_basis_type>().rmin())
                                   / ddc::discrete_space<lagrange_basis_type>().length())
                           * ddc::discrete_space<lagrange_basis_type>().length();
            }
        } else {
            if (coord_eval_interest < ddc::discrete_space<lagrange_basis_type>().rmin()) {
                return m_lower_extrap_rule(coord_eval_interest, lagrange_coef);
            }
            if (coord_eval_interest > ddc::discrete_space<lagrange_basis_type>().rmax()) {
                return m_upper_extrap_rule(coord_eval_interest, lagrange_coef);
            }
        }
        return eval_no_bc(Idx<>(), coord_eval_interest, lagrange_coef);
    }

    template <class... DerivDims, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType eval_no_bc(
            Idx<DerivDims...> const& deriv_order,
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, lagrange_domain_type, memory_space, Layout> const lagrange_coef)
            const
    {
        using deriv_dim = ddc::Deriv<continuous_dimension_type>;
        using knot_grid = std::conditional_t<
                LagrangeBasis::is_uniform(),
                UniformLagrangeKnots<LagrangeBasis>,
                NonUniformLagrangeKnots<LagrangeBasis>>;
        static_assert(
                sizeof...(DerivDims) == 0
                        || ddc::type_seq_same_v<
                                ddc::detail::TypeSeq<DerivDims...>,
                                ddc::detail::TypeSeq<deriv_dim>>,
                "The only valid dimension for deriv_order is Deriv<Dim>");

        std::array<DataType, lagrange_basis_type::degree() + 1> vals_ptr;
        Kokkos::mdspan<
                DataType,
                Kokkos::extents<std::size_t, lagrange_basis_type::degree() + 1>> const
                vals(vals_ptr.data());
        Coord<continuous_dimension_type> const coord_eval_interest(coord_eval);
        Idx<knot_grid> closest_knot = getclosest(coord_eval);
        Idx<knot_grid> first_knot
                = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        Idx<knot_grid> last_knot = ddc::discrete_space<LagrangeBasis>().break_point_domain().back();
        // The step from the closest knot to the start of the domain
        IdxStep<knot_grid> step(
                -static_cast<int>(lagrange_basis_type::degree() / 2) // back over half the domain
                - static_cast<int>(
                          lagrange_basis_type::degree()
                          % 2) // ensure coord is in the central cell if even number of points
                          * static_cast<int>(ddc::coordinate(closest_knot) > coord_eval));
        Idx<knot_grid> first_lagrange_knot;

        if constexpr (lagrange_basis_type::is_periodic()) {
            first_lagrange_knot = closest_knot
                                  + (last_knot - first_knot)
                                            * static_cast<int>((first_knot - closest_knot) > step)
                                  + step;
        } else {
            if ((first_knot - closest_knot) > step) {
                // if first_lagrange_knot would be < 0
                first_lagrange_knot = first_knot;
            } else {
                first_lagrange_knot = closest_knot + step;
                Idx<knot_grid> last_possible = last_knot - (lagrange_basis_type::degree());
                if (first_lagrange_knot > last_possible) {
                    first_lagrange_knot = last_possible;
                }
            }
        }

        static_assert(sizeof...(DerivDims) == 0);
        ddc::discrete_space<lagrange_basis_type>()
                .eval_basis(vals, coord_eval_interest, first_lagrange_knot);

        DataType y = 0.0;
        for (std::size_t i = 0; i < lagrange_basis_type::degree() + 1; ++i) {
            y += lagrange_coef(first_lagrange_knot + i) * vals[i];
        }
        return y;
    }

    /** @brief Get the knot index closest to the provided coordinate.
     */
    KOKKOS_INLINE_FUNCTION auto getclosest(Coord<continuous_dimension_type> x_interp) const
    {
        using knot_grid = std::conditional_t<
                LagrangeBasis::is_uniform(),
                UniformLagrangeKnots<LagrangeBasis>,
                NonUniformLagrangeKnots<LagrangeBasis>>;
        Idx<knot_grid> first = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        if constexpr (ddc::is_uniform_point_sampling_v<InterpolationGrid>) {
            int knot_offset = static_cast<int>(
                    Kokkos::round((x_interp - ddc::coordinate(first)) / ddc::step<knot_grid>()));
            return first + knot_offset;
        } else {
            Idx<knot_grid> last = ddc::discrete_space<LagrangeBasis>().break_point_domain().back();
            KOKKOS_ASSERT(x_interp >= ddc::discrete_space<LagrangeBasis>().rmin());
            KOKKOS_ASSERT(x_interp <= ddc::discrete_space<LagrangeBasis>().rmax());
            Idx<knot_grid> elm_cell = first + (last - first) / 2;
            while (x_interp < ddc::coordinate(elm_cell)
                   || x_interp > ddc::coordinate(elm_cell + IdxStep<knot_grid>(1))) {
                if (x_interp < ddc::coordinate(elm_cell)) {
                    last = elm_cell;
                } else {
                    first = elm_cell;
                }

                elm_cell = first + (last - first) / 2;
            }
            return elm_cell;
        }
    }
};
