// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

/**
 * @brief An ND builder that copies function values directly as interpolation coefficients.
 *
 * An ND extension of IdentityInterpolationBuilder. No computation is required because
 * the interpolation coefficients equal the function values at the grid nodes (e.g. for
 * an ND Lagrange interpolation on a tensor-product grid of nodes).
 *
 * @tparam ExecSpace           The Kokkos execution space.
 * @tparam MemorySpace         The Kokkos memory space.
 * @tparam DataType            The data type of field values and coefficients.
 * @tparam InterpolationIdxRange  The ND index range for the interpolation mesh,
 *                             of the form IdxRange<Grid1, Grid2, ...>.
 * @tparam Bases               The basis types, one per interpolation dimension, in the
 *                             same order as the grids in InterpolationIdxRange.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class InterpolationIdxRange,
        class... Bases>
class NDIdentityInterpolationBuilder;

template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class... InterpolationGrids,
        class... Bases>
class NDIdentityInterpolationBuilder<
        ExecSpace,
        MemorySpace,
        DataType,
        IdxRange<InterpolationGrids...>,
        Bases...>
{
    static_assert(sizeof...(InterpolationGrids) == sizeof...(Bases));
    static_assert(sizeof...(InterpolationGrids) > 0);

public:
    /// @brief The type of the Kokkos execution space.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space.
    using memory_space = MemorySpace;

    /// @brief The data type that the data is saved on.
    using data_type = DataType;

    /// @brief The ND index range for the interpolation mesh.
    using interpolation_idx_range_type = IdxRange<InterpolationGrids...>;

    using coeff_idx_range_type
            = IdxRange<typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>;

    /**
     * @brief Batched domain with each interpolation grid replaced by its basis domain.
     *
     * Chains ddc::replace_dim_of_t for each (InterpolationGrid_i, BasisDomain_i) pair.
     */
    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t<
                    ddc::to_type_seq_t<BatchedInterpolationIdxRange>,
                    ddc::detail::TypeSeq<InterpolationGrids...>,
                    ddc::detail::TypeSeq<
                            typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>>>;

public:
    NDIdentityInterpolationBuilder() = default;

    /**
     * @brief Compute the interpolation coefficients for a function.
     *
     * Copies vals directly to coeffs. No computation is needed because the interpolation
     * coefficients equal the function values at the grid nodes.
     *
     * @note Periodic bases are not yet supported. The copy logic for periodic dimensions
     *       (wrapping the last point back to the start) is non-trivial in ND because
     *       each combination of periodic/non-periodic dimensions must be handled
     *       separately (2^k copies for k periodic dimensions). This will be implemented
     *       when needed.
     *
     * @param[out] coeffs The coefficients of the interpolation.
     * @param[in]  vals   The values of the function on the interpolation mesh.
     */
    template <class BatchedInterpolationIdxRange>
    void operator()(
            Field<DataType,
                  batched_basis_idx_range_type<BatchedInterpolationIdxRange>,
                  memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals) const
    {
        static_assert(
                (!Bases::is_periodic() && ...),
                "NDIdentityInterpolationBuilder: periodic bases are not yet supported. "
                "For periodic dimensions the coefficient array has one extra wrap-around "
                "point per periodic dimension; implement the 2^k-copy logic and remove "
                "this assert.");
        Kokkos::deep_copy(coeffs.allocation_kokkos_view(), vals.allocation_kokkos_view());
    }

    /**
     * @brief Get the batched basis index range for a given batched interpolation domain.
     *
     * @param batched_interpolation_domain The full batched interpolation domain.
     * @return The batched basis index range.
     */
    template <class BatchedInterpolationIdxRange>
    batched_basis_idx_range_type<BatchedInterpolationIdxRange> batched_basis_idx_range(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        return batched_basis_idx_range_type<BatchedInterpolationIdxRange>(
                ddc::discrete_space<Bases>().full_domain()...,
                batched_interpolation_domain);
    }

private:
//template<class ChosenGridDim>
//copy_periodic_data(IdxRangeFull) {
//    using Dim = typename ChosenGridDim::continuous_dimension_type;
//    if constexpr (Dim::PERIODIC) {
//        using Basis = find_grid_t<Dim, ddc::detail::TypeSeq<Bases...>>;
//        using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;
//        IdxRange<basis_domain_type> bp_idx_range
//                = ddc::discrete_space<Basis>().break_point_domain().remove_last(
//                        IdxStep<basis_domain_type>(static_cast<int>(Basis::is_periodic())));
//        IdxRange<basis_domain_type> extended_domain(
//                ddc::discrete_space<Basis>().full_domain().remove_first(
//                        bp_idx_range.extents()));
//        typename BatchedInterpolationIdxRange::discrete_vector_type nrepeat(
//                extended_domain.size());
//        BatchedInterpolationIdxRange repeat_domain(get_idx_range(vals).take_first(nrepeat));
//        Kokkos::deep_copy(
//                coeffs[extended_domain].allocation_kokkos_view(),
//                vals[repeat_domain].allocation_kokkos_view());
//    }
//}
};
