#pragma once

/**
 * @brief A class to evaluate, differentiate or integrate a Lagrange function.
 *
 * A class which contains an operator () which can be used to evaluate, differentiate or integrate a Lagrange function.
 *
 * @tparam ExecSpace The Kokkos execution space on which the Lagrange evaluation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (Lagrange coefficients and evaluation) is stored.
 * @tparam BSplines The discrete dimension representing the B-Lagranges.
 * @tparam EvaluationDDim The discrete dimension on which evaluation points are defined.
 * @tparam LowerExtrapolationRule The lower extrapolation rule type.
 * @tparam UpperExtrapolationRule The upper extrapolation rule type.
 */
template <
        class ExecSpace,
        class MemorySpace,
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

    /// @brief The discrete dimension representing the B-Lagranges.
    using lagrange_basis_type = LagrangeBasis;

    /**
     * @brief The type of the whole domain representing evaluation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_evaluation_domain_type = BatchedInterpolationIdxRange;

    /// @brief The type of the 1D Lagrange domain corresponding to the dimension of interest.
    using interpolation_domain_type = IdxRange<InterpolationGrid>;

    /**
     * @brief The type of the whole Lagrange domain (cartesian product of 1D Lagrange domain
     * and batch domain) preserving the order of dimensions.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_lagrange_domain_type = BatchedInterpolationGrid;

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
                    ConstField<DataType, interpolation_domain_type, memory_space>>,
            "LowerExtrapolationRule::operator() has to be callable with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    DataType,
                    UpperExtrapolationRule,
                    Coord<continuous_dimension_type>,
                    ConstField<DataType, interpolation_domain_type, memory_space>>,
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
     * Extrapolation rules are functors used to define the behavior of the LagrangeEvaluator out of the domain where the break points of the B-Lagranges are defined.
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
     * Extrapolation rules are functors used to define the behavior of the LagrangeEvaluator out of the domain where the break points of the B-Lagranges are defined.
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
     * @brief Evaluate 1D Lagrange function (described by its Lagrange coefficients) at a given coordinate.
     *
     * The Lagrange coefficients represent a 1D Lagrange function defined on a B-Lagranges (basis Lagranges). They can be obtained via various methods, such as using a SplineBuilder.
     *
     * Remark: calling SplineBuilder then SplineEvaluator corresponds to a Lagrange interpolation.
     *
     * @param coord_eval The coordinate where the Lagrange is evaluated. Note that only the component along the dimension of interest is used.
     * @param Lagrange_coef A ChunkSpan storing the 1D Lagrange coefficients.
     *
     * @return The value of the Lagrange function at the desired coordinate.
     */
    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION DataType operator()(
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, interpolation_domain_type, memory_space, Layout> const
                    func_vals_on_grid) const
    {
        return eval(coord_eval, Lagrange_coef);
    }

    /**
     * @brief Evaluate Lagrange function (described by its Lagrange coefficients) on a mesh.
     *
     * The Lagrange coefficients represent a Lagrange function defined on a cartesian product of batch_domain and B-Lagranges
     * (basis Lagranges). They can be obtained via various methods, such as using a SplineBuilder.
     *
     * This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates
     * identified by a batch_domain_type::discrete_element_type, the evaluation is performed with the 1D set of
     * Lagrange coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * Remark: calling SplineBuilder then SplineEvaluator corresponds to a Lagrange interpolation.
     *
     * @param[out] Lagrange_eval The values of the Lagrange function at the desired coordinates. For practical reasons those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type.
     * @param[in] coords_eval The coordinates where the Lagrange is evaluated. Those are
     * stored in a ChunkSpan defined on a batched_evaluation_domain_type. Note that the coordinates of the
     * points represented by this domain are unused and irrelevant (but the points themselves (DiscreteElement) are used to select
     * the set of 1D Lagrange coefficients retained to perform the evaluation).
     * @param[in] Lagrange_coef A ChunkSpan storing the Lagrange coefficients.
     */
    template <
            class Layout1,
            class Layout2,
            class Layout3,
            class BatchedInterpolationGrid,
            class... CoordsDims>
    void operator()(
            Field<DataType, BatchedInterpolationGrid, memory_space, Layout1> const lagrange_eval,
            ConstField<Coord<CoordsDims...>, BatchedInterpolationGrid, memory_space, Layout2> const
                    coords_eval,
            ConstField<
                    DataType,
                    batched_lagrange_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout3> const lagrange_coef) const
    {
        evaluation_domain_type const evaluation_domain(lagrange_eval.domain());
        batch_domain_type<BatchedInterpolationGrid> const batch_domain(lagrange_eval.domain());

        ddc::parallel_for_each(
                "Lagrange_evaluate",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_lAMBDA(typename batch_domain_type<
                                    BatchedInterpolationGrid>::discrete_element_type const j) {
                    auto const lagrange_eval_1D = lagrange_eval[j];
                    auto const coords_eval_1D = coords_eval[j];
                    auto const lagrange_coef_1D = lagrange_coef[j];
                    for (auto const i : evaluation_domain) {
                        lagrange_eval_1D(i) = eval(coords_eval_1D(i), lagrange_coef_1D);
                    }
                });
    }

    /**
     * @brief Evaluate a spline function (described by its spline coefficients) on a mesh.
     *
     * The spline coefficients represent a spline function defined on a cartesian
     * product of batch_domain and B-splines (basis splines). They can be obtained
     * via various methods, such as using a SplineBuilder.
     *
     * This is not a multidimensional evaluation. This is a batched 1D evaluation.
     * This means that for each slice of spline_eval the evaluation is performed with
     * the 1D set of spline coefficients identified by the same batch_domain_type::discrete_element_type.
     *
     * Remark: calling SplineBuilder then SplineEvaluator corresponds to a spline
     * interpolation.
     *
     * @param[out] spline_eval The values of the spline function at the coordinates
     * of the mesh.
     * @param[in] spline_coef A ChunkSpan storing the spline coefficients.
     */
    template <class Layout1, class Layout2, class BatchedInterpolationGrid>
    void operator()(
            Field<DataType, BatchedInterpolationGrid, memory_space, Layout1> const spline_eval,
            ConstField<
                    DataType,
                    batched_spline_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout2> const spline_coef) const
    {
        evaluation_domain_type const evaluation_domain(spline_eval.domain());
        batch_domain_type<BatchedInterpolationGrid> const batch_domain(spline_eval.domain());

        ddc::parallel_for_each(
                "ddc_splines_evaluate",
                exec_space(),
                batch_domain,
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type<
                                    BatchedInterpolationGrid>::discrete_element_type const j) {
                    auto const spline_eval_1D = spline_eval[j];
                    auto const spline_coef_1D = spline_coef[j];
                    for (auto const i : evaluation_domain) {
                        Coord<continuous_dimension_type> coord_eval_1D = ddc::coordinate(i);
                        spline_eval_1D(i) = eval(coord_eval_1D, spline_coef_1D);
                    }
                });
    }

private:
    template <class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType
    eval(Coord<CoordsDims...> const& coord_eval,
         ConstField<DataType, spline_domain_type, memory_space, Layout> const spline_coef) const
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
                return m_lower_extrap_rule(coord_eval_interest, spline_coef);
            }
            if (coord_eval_interest > ddc::discrete_space<lagrange_basis_type>().rmax()) {
                return m_upper_extrap_rule(coord_eval_interest, spline_coef);
            }
        }
        return eval_no_bc(Idx<>(), coord_eval_interest, spline_coef);
    }

    template <class... DerivDims, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType eval_no_bc(
            Idx<DerivDims...> const& deriv_order,
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, spline_domain_type, Layout, memory_space> const spline_coef) const
    {
        using deriv_dim = ddc::Deriv<continuous_dimension_type>;
        using knot_grid = std::conditional_t<
                LagrangeBasis::is_uniform(),
                UniformLagrangeKnots<LagrangeBasis>,
                NonUniformLagrangeKnots<LagrangeBasis>>;
        static_assert(
                sizeof...(DerivDims) == 0
                        || type_seq_same_v<
                                detail::TypeSeq<DerivDims...>,
                                detail::TypeSeq<deriv_dim>>,
                "The only valid dimension for deriv_order is Deriv<Dim>");

        std::array<DataType, lagrange_basis_type::degree() + 1> vals_ptr;
        Kokkos::mdspan<
                DataType,
                Kokkos::extents<std::size_t, lagrange_basis_type::degree() + 1>> const
                vals(vals_ptr.data());
        Coord<continuous_dimension_type> const coord_eval_interest(coord_eval);
        Idx<knot_grid> closest_knot = getclosest(x);
        Idx<knot_grid> first_knot
                = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        Idx<knot_grid> last_knot
                = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        IdxStep<knot_grid> step(
                -lagrange_basis_type::degree() / 2
                - (lagrange_basis_type::degree() % 2)
                          * static_cast<int>(ddc::coordinate(closest_knot) > x));
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
                Idx<knot_grid> last_possible = last_knot - (lagrange_basis_type::degree() + 1);
                if (first_lagrange_knot > last_possible) {
                    first_lagrange_knot = last_possible;
                }
            }
        }

        static_assert(sizeof...(DerivDims) == 0);
        ddc::discrete_space<lagrange_basis_type>()
                .eval_basis(vals, first_lagrange_knot, coord_eval_interest);

        DataType y = 0.0;
        for (std::size_t i = 0; i < lagrange_basis_type::degree() + 1; ++i) {
            y += spline_coef(Idx<lagrange_basis_type>(jmin + i)) * vals[i];
        }
        return y;
    }

    KOKKOS_INLINE_FUNCTION auto getclosest(coord_type x) const
    {
        using knot_grid = std::conditional_t<
                LagrangeBasis::is_uniform(),
                UniformLagrangeKnots<LagrangeBasis>,
                NonUniformLagrangeKnots<LagrangeBasis>>;
        Idx<knot_grid> first = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        if constexpr (ddc::is_uniform_point_sampling_v<GridInterp>) {
            int knot_offset = static_cast<int>(Kokkos::round(
                    (x - ddc::coordinate(first))
                    / ddc::step<decltype(ddc::discrete_space<LagrangeBasis>())::knot_grid>()));
            return first + knot_offset;
        } else {
            Idx<knot_grid> last = ddc::discrete_space<LagrangeBasis>().break_point_domain().back();
            KOKKOS_ASSERT(x_interp >= ddc::discrete_space<LagrangeBasis>().rmin());
            KOKKOS_ASSERT(x_interp <= ddc::discrete_space<LagrangeBasis>().rmax());
            IdxInterp elm_cell = begin + (end - begin) / 2;
            while (x_interp < ddc::coordinate(elm_cell)
                   || x_interp > ddc::coordinate(elm_cell + IdxStepInterp(1))) {
                if (x_interp < ddc::coordinate(elm_cell)) {
                    end = elm_cell;
                } else {
                    begin = elm_cell;
                }

                elm_cell = begin + (end - begin) / 2;
            }
            return elm_cell;
        }
    }
};
