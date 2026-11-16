

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
#include "tensor_common.hpp"
#include "vector_field.hpp"

namespace detail {

template <class ElementDataType, class = void>
struct QuadratureDataType
{
    static_assert(std::is_floating_point_v<ElementDataType>);
    using data_type = ElementDataType;
};

template <class ElementDataType>
struct QuadratureDataType<ElementDataType, std::enable_if_t<is_tensor_type_v<ElementDataType>>>
{
    static_assert(ElementDataType::rank() == 1);
    using data_type = typename ElementDataType::element_type;
};

} // namespace detail

template <
        class IdxRangeQuadrature,
        class IdxRangeTotal = IdxRangeQuadrature,
        class ElementDataType = double,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
class Quadrature
{
    using DataType = typename detail::QuadratureDataType<ElementDataType>::data_type;

    static_assert(std::is_floating_point_v<DataType>);

private:
    using IdxQuadrature = typename IdxRangeQuadrature::discrete_element_type;
    using IdxStepQuadrature = typename IdxRangeQuadrature::discrete_vector_type;

    using QuadConstField = ConstField<DataType, IdxRangeQuadrature, MemorySpace>;

    QuadConstField m_coefficients;

public:
    explicit Quadrature(QuadConstField coeffs) : m_coefficients(coeffs) {}

    template <class ExecutionSpace, class IntegratorFunction>
    ElementDataType operator()(ExecutionSpace exec_space, IntegratorFunction integrated_function)
            const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        static_assert(
                std::is_invocable_v<IntegratorFunction, IdxQuadrature>,
                "The object passed to Quadrature::operator() is not defined on the quadrature "
                "idx_range.");

        QuadConstField const coeff_proxy = m_coefficients;

        if constexpr (ddc::is_chunk_v<IntegratorFunction>) {
            // Sanity check ensure that the front and back are in the domain
            // If this fails data may be distributed across MPI ranks
            assert(get_idx_range(integrated_function).contains(get_idx_range(coeff_proxy).front()));
            assert(get_idx_range(integrated_function).contains(get_idx_range(coeff_proxy).back()));
        }

        // This fence helps avoid a CPU seg fault. See #290 for more details
        exec_space.fence();
        // This condition is necessary to execute in serial, even in a device activated build.
        // Without it a seg fault appears
        if constexpr (std::is_same_v<ExecutionSpace, Kokkos::DefaultHostExecutionSpace>) {
            return ddc::host_transform_reduce(
                    get_idx_range(coeff_proxy),
                    ElementDataType(0.0),
                    ddc::reducer::sum<ElementDataType>(),
                    KOKKOS_LAMBDA(IdxQuadrature const ix) {
                        return coeff_proxy(ix) * integrated_function(ix);
                    });
        } else {
            const std::source_location location = std::source_location::current();
            return ddc::parallel_transform_reduce(
                    location.function_name(),
                    exec_space,
                    get_idx_range(coeff_proxy),
                    ElementDataType(0.0),
                    ddc::reducer::sum<ElementDataType>(),
                    KOKKOS_LAMBDA(IdxQuadrature const ix) {
                        return coeff_proxy(ix) * integrated_function(ix);
                    });
        }
    }

    template <class ExecutionSpace, class ResultField, class IntegratorFunction>
    void operator()(
            ExecutionSpace exec_space,
            ResultField const result,
            IntegratorFunction integrated_function) const
            requires(std::is_floating_point_v<
                     ElementDataType>) //TOD: Remove requires when #322 is fixed
    {
        using IdxRangeBatch = ResultField::discrete_domain_type;
        static_assert(std::is_same_v<typename ResultField::element_type, ElementDataType>);
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        using ExpectedBatchDims = ddc::type_seq_remove_t<
                ddc::to_type_seq_t<IdxRangeTotal>,
                ddc::to_type_seq_t<IdxRangeQuadrature>>;
        static_assert(
                ddc::type_seq_same_v<ddc::to_type_seq_t<IdxRangeBatch>, ExpectedBatchDims>,
                "The batch idx_range deduced from the type of result does not match the class "
                "template parameters.");

        if constexpr (ddc::is_chunk_v<IntegratorFunction>) {
            // Sanity check ensure that the front and back are in the domain
            // If this fails data may be distributed across MPI ranks
            assert(IdxRangeQuadrature(get_idx_range(integrated_function))
                           .contains(get_idx_range(m_coefficients).front()));
            assert(IdxRangeQuadrature(get_idx_range(integrated_function))
                           .contains(get_idx_range(m_coefficients).back()));
        }

        // Get useful index types
        using IdxTotal = typename IdxRangeTotal::discrete_element_type;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        static_assert(
                std::is_invocable_v<IntegratorFunction, IdxTotal>,
                "The object passed to Quadrature::operator() is not defined on the total "
                "idx_range.");

        // Get index ranges
        IdxRangeQuadrature quad_idx_range(get_idx_range(m_coefficients));
        IdxRangeBatch batch_idx_range(get_idx_range(result));

        QuadConstField const coeff_proxy = m_coefficients;
        // Loop over batch dimensions
        ddc::parallel_for_each_team(
                exec_space,
                batch_idx_range,
                KOKKOS_LAMBDA(
                        ddc::TeamPolicyMember<ExecutionSpace> const& team,
                        IdxBatch const ib) {
                    // Sum over quadrature dimensions
                    ElementDataType teamSum(0);
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, quad_idx_range.size()),
                            [&](int const& thread_index, ElementDataType& sum) {
                                IdxQuadrature iq = ddcHelper::
                                        get_idx_from_linear_index(quad_idx_range, thread_index);
                                IdxTotal it(ib, iq);
                                sum += coeff_proxy(iq) * integrated_function(it);
                            },
                            teamSum);
                    result(ib) = teamSum;
                });
    }
};

// Auto-deduce template parameters
template <class DataType, class IdxRangeQuadrature, class MemorySpace>
Quadrature(ConstField<DataType, IdxRangeQuadrature, MemorySpace>)
        -> Quadrature<IdxRangeQuadrature, IdxRangeQuadrature, DataType, MemorySpace>;

namespace detail {
template <
        class NewMemorySpace,
        class IdxRangeQuadrature,
        class IdxRangeTotal,
        class ElementDataType,
        class MemorySpace>
struct OnMemorySpace<
        NewMemorySpace,
        Quadrature<IdxRangeQuadrature, IdxRangeTotal, ElementDataType, MemorySpace>>
{
    using type = Quadrature<IdxRangeQuadrature, IdxRangeTotal, ElementDataType, NewMemorySpace>;
};
} // namespace detail
```


