

# File lagrange\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**lagrange\_evaluator.hpp**](lagrange__evaluator_8hpp.md)

[Go to the documentation of this file](lagrange__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "i_interpolation_evaluator.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"

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
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using data_type = DataType;

    using continuous_dimension_type = typename LagrangeBasis::continuous_dimension_type;

    using lagrange_basis_type = LagrangeBasis;

    using coeff_grid_type =
            typename LagrangeBasis::template Impl<LagrangeBasis, MemorySpace>::knot_grid;

    using evaluation_idx_range_type = IdxRange<InterpolationGrid>;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_evaluation_idx_range_type = BatchedInterpolationIdxRange;

    using coeff_idx_range_type = IdxRange<coeff_grid_type>;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batch_idx_range_type
            = ddc::remove_dims_of_t<BatchedInterpolationIdxRange, InterpolationGrid>;

    template <
            class BatchedInterpolationIdxRange,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange>>>
    using batched_coeff_idx_range_type = ddc::
            replace_dim_of_t<BatchedInterpolationIdxRange, InterpolationGrid, coeff_grid_type>;

    using lower_extrapolation_rule_type = LowerExtrapolationRule;

    using upper_extrapolation_rule_type = UpperExtrapolationRule;


private:
    template <class, class...>
    friend class NDLagrangeEvaluator;

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
                    ConstField<DataType, coeff_idx_range_type, memory_space>>,
            "LowerExtrapolationRule::operator() has to be callable with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    DataType,
                    UpperExtrapolationRule,
                    Coord<continuous_dimension_type>,
                    ConstField<DataType, coeff_idx_range_type, memory_space>>,
            "UpperExtrapolationRule::operator() has to be callable with usual arguments.");

    explicit LagrangeEvaluator(
            LowerExtrapolationRule const& lower_extrap_rule,
            UpperExtrapolationRule const& upper_extrap_rule)
        : m_lower_extrap_rule(lower_extrap_rule)
        , m_upper_extrap_rule(upper_extrap_rule)
    {
    }

    LagrangeEvaluator(LagrangeEvaluator const& x) = default;

    LagrangeEvaluator(LagrangeEvaluator&& x) = default;

    ~LagrangeEvaluator() = default;

    LagrangeEvaluator& operator=(LagrangeEvaluator const& x) = default;

    LagrangeEvaluator& operator=(LagrangeEvaluator&& x) = default;

    lower_extrapolation_rule_type lower_extrapolation_rule() const
    {
        return m_lower_extrap_rule;
    }

    upper_extrapolation_rule_type upper_extrapolation_rule() const
    {
        return m_upper_extrap_rule;
    }

    template <class Layout, class... CoordsDims, class BatchedLagrangeIdxRange>
    KOKKOS_FUNCTION DataType operator()(
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, BatchedLagrangeIdxRange, memory_space, Layout> const lagrange_coef)
            const
    {
        return eval(coord_eval, lagrange_coef);
    }

    template <
            class Layout1,
            class Layout2,
            class Layout3,
            class IdxRangeBatchedInterpolation,
            class... CoordsDims>
    void operator()(
            Field<DataType, IdxRangeBatchedInterpolation, memory_space, Layout1> const
                    lagrange_eval,
            ConstField<
                    Coord<CoordsDims...>,
                    IdxRangeBatchedInterpolation,
                    memory_space,
                    Layout2> const coords_eval,
            ConstField<
                    DataType,
                    batched_coeff_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout3> const lagrange_coef) const
    {
        evaluation_idx_range_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_idx_range_type<IdxRangeBatchedInterpolation> const batch_idx_range(
                get_idx_range(lagrange_eval));

        using IdxBatchInterpolation =
                typename batch_idx_range_type<IdxRangeBatchedInterpolation>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(IdxBatchInterpolation const j) {
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        lagrange_eval(j, i) = eval(coords_eval(j, i), lagrange_coef[j]);
                    }
                });
    }

    template <class Layout1, class Layout2, class BatchedInterpolationIdxRange>
    void operator()(
            Field<DataType, BatchedInterpolationIdxRange, memory_space, Layout1> const
                    lagrange_eval,
            ConstField<
                    DataType,
                    batched_coeff_idx_range_type<BatchedInterpolationIdxRange>,
                    memory_space,
                    Layout2> const lagrange_coef) const
    {
        evaluation_idx_range_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_idx_range_type<BatchedInterpolationIdxRange> const batch_idx_range(
                get_idx_range(lagrange_eval));

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(typename batch_idx_range_type<
                                    BatchedInterpolationIdxRange>::discrete_element_type const j) {
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        Coord<continuous_dimension_type> coord_eval_1D = ddc::coordinate(i);
                        lagrange_eval(j, i) = eval(coord_eval_1D, lagrange_coef[j]);
                    }
                });
    }

    template <class IdxDerivDims, class Layout, class BatchedLagrangeIdxRange, class... CoordsDims>
    KOKKOS_FUNCTION double deriv(
            IdxDerivDims const& deriv_order,
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, BatchedLagrangeIdxRange, memory_space, Layout> const lagrange_coef)
            const
    {
        return eval_no_bc(deriv_order, coord_eval, lagrange_coef);
    }

    template <
            class IdxDerivDims,
            class Layout1,
            class Layout2,
            class Layout3,
            class IdxRangeBatchedInterpolation,
            class... CoordsDims>
    void deriv(
            IdxDerivDims const& deriv_order,
            Field<DataType, IdxRangeBatchedInterpolation, memory_space, Layout1> const
                    lagrange_eval,
            ConstField<
                    Coord<CoordsDims...>,
                    IdxRangeBatchedInterpolation,
                    memory_space,
                    Layout2> const coords_eval,
            ConstField<
                    DataType,
                    batched_coeff_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout3> const lagrange_coef) const
    {
        evaluation_idx_range_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_idx_range_type<IdxRangeBatchedInterpolation> const batch_idx_range(
                get_idx_range(lagrange_eval));

        using IdxBatchInterpolation =
                typename batch_idx_range_type<IdxRangeBatchedInterpolation>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(IdxBatchInterpolation const j) {
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        lagrange_eval(j, i)
                                = eval_no_bc(deriv_order, coords_eval(j, i), lagrange_coef[j]);
                    }
                });
    }

    template <class IdxDerivDims, class Layout1, class Layout2, class IdxRangeBatchedInterpolation>
    void deriv(
            IdxDerivDims const& deriv_order,
            Field<DataType, IdxRangeBatchedInterpolation, memory_space, Layout1> const
                    lagrange_eval,
            ConstField<
                    DataType,
                    batched_coeff_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout2> const lagrange_coef) const
    {
        evaluation_idx_range_type const evaluation_idx_range(get_idx_range(lagrange_eval));
        batch_idx_range_type<IdxRangeBatchedInterpolation> const batch_idx_range(
                get_idx_range(lagrange_eval));

        using IdxBatchInterpolation =
                typename batch_idx_range_type<IdxRangeBatchedInterpolation>::discrete_element_type;

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space(),
                batch_idx_range,
                KOKKOS_CLASS_LAMBDA(IdxBatchInterpolation const j) {
                    for (Idx<InterpolationGrid> const i : evaluation_idx_range) {
                        lagrange_eval(j, i)
                                = eval_no_bc(deriv_order, ddc::coordinate(i), lagrange_coef[j]);
                    }
                });
    }

private:
    template <class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType
    eval(Coord<CoordsDims...> const& coord_eval,
         ConstField<DataType, coeff_idx_range_type, memory_space, Layout> const lagrange_coef) const
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

    KOKKOS_INLINE_FUNCTION Idx<coeff_grid_type> find_stencil(
            Coord<continuous_dimension_type>& coord) const
    {
        Idx<coeff_grid_type> const closest_knot = getclosest(coord);
        Idx<coeff_grid_type> const first_knot
                = ddc::discrete_space<lagrange_basis_type>().break_point_domain().front();
        Idx<coeff_grid_type> const last_knot
                = ddc::discrete_space<lagrange_basis_type>().break_point_domain().back();
        // The step from the closest knot to the start of the stencil
        IdxStep<coeff_grid_type> const step(
                -static_cast<int>(lagrange_basis_type::degree() / 2)
                - static_cast<int>(lagrange_basis_type::degree() % 2)
                          * static_cast<int>(ddc::coordinate(closest_knot) > coord));
        Idx<coeff_grid_type> first_lagrange_knot;
        if constexpr (lagrange_basis_type::is_periodic()) {
            first_lagrange_knot = closest_knot
                                  + (last_knot - first_knot)
                                            * static_cast<int>((first_knot - closest_knot) > step)
                                  + step;
            if (first_lagrange_knot > closest_knot) {
                coord += ddc::discrete_space<lagrange_basis_type>().length();
            }
        } else {
            if ((first_knot - closest_knot) > step) {
                first_lagrange_knot = first_knot;
            } else {
                first_lagrange_knot = closest_knot + step;
                Idx<coeff_grid_type> const last_possible
                        = last_knot - lagrange_basis_type::degree();
                if (first_lagrange_knot > last_possible) {
                    first_lagrange_knot = last_possible;
                }
            }
        }
        return first_lagrange_knot;
    }

    template <class... DerivDims, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION DataType eval_no_bc(
            Idx<DerivDims...> const& deriv_order,
            Coord<CoordsDims...> const& coord_eval,
            ConstField<DataType, coeff_idx_range_type, memory_space, Layout> const lagrange_coef)
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

        constexpr std::size_t n_basis = lagrange_basis_type::degree() + 1;
        std::array<DataType, n_basis> vals_ptr;
        Kokkos::mdspan<DataType, Kokkos::extents<std::size_t, n_basis>> const vals(vals_ptr.data());

        Coord<continuous_dimension_type> coord_interest(coord_eval);
        Idx<knot_grid> first_lagrange_knot = find_stencil(coord_interest);
        ddc::discrete_space<lagrange_basis_type>()
                .eval_basis(vals, coord_interest, first_lagrange_knot);
        if constexpr (sizeof...(DerivDims) == 0) {
            ddc::discrete_space<lagrange_basis_type>()
                    .eval_basis(vals, coord_interest, first_lagrange_knot);
        } else {
            unsigned long const order = deriv_order.uid();
            KOKKOS_ASSERT(order <= lagrange_basis_type::degree());

            std::array<DataType, n_basis * n_basis> derivs_ptr;
            Kokkos::mdspan<
                    DataType,
                    Kokkos::extents<std::size_t, n_basis, Kokkos::dynamic_extent>> const
                    derivs(derivs_ptr.data(), order + 1);

            ddc::discrete_space<lagrange_basis_type>()
                    .eval_basis_and_n_derivs(derivs, coord_interest, first_lagrange_knot, order);

            for (std::size_t i = 0; i < n_basis; ++i) {
                vals(i) = derivs(i, order);
            }
        }

        DataType y = 0.0;
        for (std::size_t i = 0; i < n_basis; ++i) {
            y += lagrange_coef(first_lagrange_knot + i) * vals[i];
        }
        return y;
    }

    KOKKOS_INLINE_FUNCTION auto getclosest(Coord<continuous_dimension_type> x_interp) const
    {
        using knot_grid = std::conditional_t<
                LagrangeBasis::is_uniform(),
                UniformLagrangeKnots<LagrangeBasis>,
                NonUniformLagrangeKnots<LagrangeBasis>>;
        Idx<knot_grid> first = ddc::discrete_space<LagrangeBasis>().break_point_domain().front();
        if constexpr (ddc::is_uniform_point_sampling_v<knot_grid>) {
            if constexpr (lagrange_basis_type::degree() % 2 == 0) {
                int knot_offset = static_cast<int>(Kokkos::round(
                        (x_interp - ddc::coordinate(first)) / ddc::step<knot_grid>()));
                return first + knot_offset;
            } else {
                int knot_offset = static_cast<int>(
                        (x_interp - ddc::coordinate(first)) / ddc::step<knot_grid>());
                return first + knot_offset;
            }
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
            DataType cell_len
                    = ddc::coordinate(elm_cell + IdxStep<knot_grid>(1)) - ddc::coordinate(elm_cell);
            return elm_cell
                   + IdxStep<knot_grid>(
                           ((lagrange_basis_type::degree() + 1) % 2)
                           * static_cast<int>(x_interp - ddc::coordinate(elm_cell) > cell_len / 2));
        }
    }
};
```


