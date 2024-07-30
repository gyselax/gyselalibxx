// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

/**
 * @brief A class providing an operator for integrating functions defined on a discrete index range.
 *
 * @tparam QuadratureIdxRange The index range over which the function is integrated.
 * @tparam TotalIdxRange The index range of the chunk which can be passed to the operator(). This is the
 *                      QuadratureIdxRange combined with any batch dimensions. If there are no
 *                      batch dimensions then this argument does not need to be provided as by
 *                      default it is equal to the QuadratureIdxRange.
 * @tparam MemorySpace The memory space (cpu/gpu) where the quadrature coefficients are saved.
 */
template <
        class QuadratureIdxRange,
        class TotalIdxRange = QuadratureIdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
class Quadrature
{
private:
    /// The tyoe of an element of an index of the quadrature coefficients.
    using QuadratureIdx = typename QuadratureIdxRange::discrete_element_type;

    using QuadConstField
            = DConstField<QuadratureIdxRange, std::experimental::layout_right, MemorySpace>;

    QuadConstField m_coefficients;

public:
    /**
     * @brief Create a Quadrature object.
     * @param coeffs
     * 	      The coefficients of the quadrature.
     */
    explicit Quadrature(QuadConstField coeffs) : m_coefficients(coeffs) {}

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete index range.
     *
     * @param[in] exec_space
     *        The space on which the function is executed (CPU/GPU).
     * @param[in] integrated_function
     *        A function taking an index of a position in the index range over which the quadrature is
     *        calculated and returning the value of the function to be integrated at that point.
     *        It should be noted that a ChunkSpan fulfils these criteria and can be passed as the function to be integrated.
     *        If the exec_space is a GPU the function that is passed must be accessible from GPU.
     *
     * @returns The integral of the function over the index range.
     */
    template <class ExecutionSpace, class IntegratorFunction>
    double operator()(ExecutionSpace exec_space, IntegratorFunction integrated_function) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        static_assert(
                std::is_invocable_v<IntegratorFunction, QuadratureIdx>,
                "The object passed to Quadrature::operator() is not defined on the quadrature "
                "idx_range.");

        ddc::ChunkSpan const coeff_proxy = m_coefficients;

        // This fence helps avoid a CPU seg fault. See #290 for more details
        exec_space.fence();
        return ddc::parallel_transform_reduce(
                exec_space,
                get_idx_range(coeff_proxy),
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(QuadratureIdx const ix) {
                    return coeff_proxy(ix) * integrated_function(ix);
                });
    }

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete index range
     * by cycling over batch dimensions.
     *
     * @param[in] exec_space
     *        The space on which the function is executed (CPU/GPU).
     * @param[out] result
     *        The result of the quadrature calculation.
     * @param[in] integrated_function
     *        A function taking an index of a position in the index range over which the quadrature is
     *        calculated (including the batch index range) and returning the value of the function to
     *        be integrated at that point.
     *        Please note that a ChunkSpan fulfills the described criteria.
     *        If the exec_space is a GPU the function that is passed must be accessible from GPU.
     */
    template <class ExecutionSpace, class BatchIdxRange, class IntegratorFunction>
    void operator()(
            ExecutionSpace exec_space,
            Field<double, BatchIdxRange, std::experimental::layout_right, MemorySpace> const result,
            IntegratorFunction integrated_function) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        static_assert(
                std::is_same_v<ExecutionSpace, Kokkos::DefaultExecutionSpace>,
                "Kokkos::TeamPolicy only works with the default execution space. Please use "
                "DefaultExecutionSpace to call this batched operator.");
        using ExpectedBatchDims = ddc::type_seq_remove_t<
                ddc::to_type_seq_t<TotalIdxRange>,
                ddc::to_type_seq_t<QuadratureIdxRange>>;
        static_assert(
                ddc::type_seq_same_v<ddc::to_type_seq_t<BatchIdxRange>, ExpectedBatchDims>,
                "The batch idx_range deduced from the type of result does not match the class "
                "template parameters.");

        // Get useful index types
        using TotalIdx = typename TotalIdxRange::discrete_element_type;
        using BatchIdx = typename BatchIdxRange::discrete_element_type;

        static_assert(
                std::is_invocable_v<IntegratorFunction, TotalIdx>,
                "The object passed to Quadrature::operator() is not defined on the total "
                "idx_range.");

        // Get index ranges
        QuadratureIdxRange quad_idx_range(get_idx_range(m_coefficients));
        BatchIdxRange batch_idx_range(get_idx_range(result));

        ddc::ChunkSpan const coeff_proxy = m_coefficients;
        // Loop over batch dimensions
        Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(exec_space, batch_idx_range.size(), Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int idx = team.league_rank();
                    BatchIdx ib = to_discrete_element(idx, batch_idx_range);

                    // Sum over quadrature dimensions
                    double teamSum = 0;
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, quad_idx_range.size()),
                            [&](int const& thread_index, double& sum) {
                                QuadratureIdx iq
                                        = to_discrete_element(thread_index, quad_idx_range);
                                TotalIdx it(ib, iq);
                                sum += coeff_proxy(iq) * integrated_function(it);
                            },
                            teamSum);
                    result(ib) = teamSum;
                });
    }

private:
    /**
     * A function which converts an integer into a DiscreteElement found in an index range
     * starting from the front. This is useful for iterating over an index range using Kokkos
     * loops.
     *
     * @param[in] idx The index of the requested element.
     * @param[in] dom The index range being iterated over.
     *
     * @return The vector displacement from the front of the index range
     */
    template <class HeadDim, class... Grid1D>
    KOKKOS_FUNCTION static Idx<HeadDim, Grid1D...> to_discrete_element(
            int idx,
            IdxRange<HeadDim, Grid1D...> dom)
    {
        IdxRange<Grid1D...> subidx_range(dom);
        Idx<HeadDim> head_idx(ddc::select<HeadDim>(dom).front() + idx / subidx_range.size());
        if constexpr (sizeof...(Grid1D) == 0) {
            return head_idx;
        } else {
            Idx<Grid1D...> tail_idx = to_discrete_element(idx % subidx_range.size(), subidx_range);
            return Idx<HeadDim, Grid1D...>(head_idx, tail_idx);
        }
    }
};

namespace detail {
template <class NewMemorySpace, class QuadratureIdxRange, class TotalIdxRange, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, Quadrature<QuadratureIdxRange, TotalIdxRange, MemorySpace>>
{
    using type = Quadrature<QuadratureIdxRange, TotalIdxRange, NewMemorySpace>;
};
} // namespace detail
