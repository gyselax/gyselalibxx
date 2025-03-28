

# File quadrature.hpp

[**File List**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**quadrature.hpp**](quadrature_8hpp.md)

[Go to the documentation of this file](quadrature_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

template <
        class IdxRangeQuadrature,
        class IdxRangeTotal = IdxRangeQuadrature,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
class Quadrature
{
private:
    using IdxQuadrature = typename IdxRangeQuadrature::discrete_element_type;

    using QuadConstField = DConstField<IdxRangeQuadrature, MemorySpace>;

    QuadConstField m_coefficients;

public:
    explicit Quadrature(QuadConstField coeffs) : m_coefficients(coeffs) {}

    template <class ExecutionSpace, class IntegratorFunction>
    double operator()(ExecutionSpace exec_space, IntegratorFunction integrated_function) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        static_assert(
                std::is_invocable_v<IntegratorFunction, IdxQuadrature>,
                "The object passed to Quadrature::operator() is not defined on the quadrature "
                "idx_range.");

        QuadConstField const coeff_proxy = m_coefficients;

        // This fence helps avoid a CPU seg fault. See #290 for more details
        exec_space.fence();
        // This condition is necessary to execute in serial, even in a device activated build.
        // Without it a seg fault appears
        if constexpr (std::is_same_v<ExecutionSpace, Kokkos::DefaultHostExecutionSpace>) {
            return ddc::transform_reduce(
                    get_idx_range(coeff_proxy),
                    0.0,
                    ddc::reducer::sum<double>(),
                    KOKKOS_LAMBDA(IdxQuadrature const ix) {
                        return coeff_proxy(ix) * integrated_function(ix);
                    });
        } else {
            return ddc::parallel_transform_reduce(
                    exec_space,
                    get_idx_range(coeff_proxy),
                    0.0,
                    ddc::reducer::sum<double>(),
                    KOKKOS_LAMBDA(IdxQuadrature const ix) {
                        return coeff_proxy(ix) * integrated_function(ix);
                    });
        }
    }

    template <class ExecutionSpace, class BatchIdxRange, class IntegratorFunction>
    void operator()(
            ExecutionSpace exec_space,
            Field<double, BatchIdxRange, MemorySpace> const result,
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
                ddc::to_type_seq_t<IdxRangeTotal>,
                ddc::to_type_seq_t<IdxRangeQuadrature>>;
        static_assert(
                ddc::type_seq_same_v<ddc::to_type_seq_t<BatchIdxRange>, ExpectedBatchDims>,
                "The batch idx_range deduced from the type of result does not match the class "
                "template parameters.");

        // Get useful index types
        using IdxTotal = typename IdxRangeTotal::discrete_element_type;
        using IdxBatch = typename BatchIdxRange::discrete_element_type;

        static_assert(
                std::is_invocable_v<IntegratorFunction, IdxTotal>,
                "The object passed to Quadrature::operator() is not defined on the total "
                "idx_range.");

        // Get index ranges
        IdxRangeQuadrature quad_idx_range(get_idx_range(m_coefficients));
        BatchIdxRange batch_idx_range(get_idx_range(result));

        QuadConstField const coeff_proxy = m_coefficients;
        // Loop over batch dimensions
        Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(exec_space, batch_idx_range.size(), Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int idx = team.league_rank();
                    IdxBatch ib = to_discrete_element(idx, batch_idx_range);

                    // Sum over quadrature dimensions
                    double teamSum = 0;
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, quad_idx_range.size()),
                            [&](int const& thread_index, double& sum) {
                                IdxQuadrature iq
                                        = to_discrete_element(thread_index, quad_idx_range);
                                IdxTotal it(ib, iq);
                                sum += coeff_proxy(iq) * integrated_function(it);
                            },
                            teamSum);
                    result(ib) = teamSum;
                });
    }

private:
    template <class HeadDim, class... Grid1D>
    KOKKOS_FUNCTION static Idx<HeadDim, Grid1D...> to_discrete_element(
            int idx,
            IdxRange<HeadDim, Grid1D...> idx_range)
    {
        IdxRange<Grid1D...> subidx_range(idx_range);
        Idx<HeadDim> head_idx(ddc::select<HeadDim>(idx_range).front() + idx / subidx_range.size());
        if constexpr (sizeof...(Grid1D) == 0) {
            return head_idx;
        } else {
            Idx<Grid1D...> tail_idx = to_discrete_element(idx % subidx_range.size(), subidx_range);
            return Idx<HeadDim, Grid1D...>(head_idx, tail_idx);
        }
    }
};

namespace detail {
template <class NewMemorySpace, class IdxRangeQuadrature, class IdxRangeTotal, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, Quadrature<IdxRangeQuadrature, IdxRangeTotal, MemorySpace>>
{
    using type = Quadrature<IdxRangeQuadrature, IdxRangeTotal, NewMemorySpace>;
};
} // namespace detail
```


