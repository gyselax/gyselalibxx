// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "ddc_helper.hpp"

/**
 * @brief A class providing an operator for integrating functions defined on a discrete domain.
 *
 * @tparam QuadratureDomain The domain over which the function is integrated.
 * @tparam TotalDomain The domain of the chunk which can be passed to the operator(). This is the
 *                      QuadratureDomain combined with any batch dimensions. If there are no
 *                      batch dimensions then this argument does not need to be provided as by
 *                      default it is equal to the QuadratureDomain.
 * @tparam MemorySpace The memory space (cpu/gpu) where the quadrature coefficients are saved.
 */
template <
        class QuadratureDomain,
        class TotalDomain = QuadratureDomain,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
class Quadrature
{
private:
    /// The tyoe of an element of an index of the quadrature coefficients.
    using QuadratureIndex = typename QuadratureDomain::discrete_element_type;

    using QuadChunkView = ddc::
            ChunkView<double, QuadratureDomain, std::experimental::layout_right, MemorySpace>;

    QuadChunkView m_coefficients;

public:
    /**
     * @brief Create a Quadrature object.
     * @param coeffs
     * 	      The coefficients of the quadrature.
     */
    explicit Quadrature(QuadChunkView coeffs) : m_coefficients(coeffs) {}

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete domain.
     *
     * @param[in] exec_space
     *        The space on which the function is executed (CPU/GPU).
     * @param[in] integrated_function
     *        A function taking an index of a position in the domain over which the quadrature is
     *        calculated and returning the value of the function to be integrated at that point.
     *        It should be noted that a ChunkSpan fulfils these criteria and can be passed as the function to be integrated.
     *        If the exec_space is a GPU the function that is passed must be accessible from GPU.
     *
     * @returns The integral of the function over the domain.
     */
    template <class ExecutionSpace, class IntegratorFunction>
    double operator()(ExecutionSpace exec_space, IntegratorFunction integrated_function) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible,
                "Execution space is not compatible with memory space where coefficients are found");
        static_assert(
                std::is_invocable_v<IntegratorFunction, QuadratureIndex>,
                "The object passed to Quadrature::operator() is not defined on the quadrature "
                "domain.");

        ddc::ChunkSpan const coeff_proxy = m_coefficients;

        // This fence helps avoid a CPU seg fault. See #290 for more details
        exec_space.fence();
        return ddc::parallel_transform_reduce(
                exec_space,
                coeff_proxy.domain(),
                0.0,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(QuadratureIndex const ix) {
                    return coeff_proxy(ix) * integrated_function(ix);
                });
    }

    /**
     * @brief An operator for calculating the integral of a function defined on a discrete domain
     * by cycling over batch dimensions.
     *
     * @param[in] exec_space
     *        The space on which the function is executed (CPU/GPU).
     * @param[out] result
     *        The result of the quadrature calculation.
     * @param[in] integrated_function
     *        A function taking an index of a position in the domain over which the quadrature is
     *        calculated (including the batch domain) and returning the value of the function to
     *        be integrated at that point.
     *        Please note that a ChunkSpan fulfills the described criteria.
     *        If the exec_space is a GPU the function that is passed must be accessible from GPU.
     */
    template <class ExecutionSpace, class BatchDomain, class IntegratorFunction>
    void operator()(
            ExecutionSpace exec_space,
            ddc::ChunkSpan<double, BatchDomain, std::experimental::layout_right, MemorySpace> const
                    result,
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
                ddc::to_type_seq_t<TotalDomain>,
                ddc::to_type_seq_t<QuadratureDomain>>;
        static_assert(
                ddc::type_seq_same_v<ddc::to_type_seq_t<BatchDomain>, ExpectedBatchDims>,
                "The batch domain deduced from the type of result does not match the class "
                "template parameters.");

        // Get useful index types
        using TotalIndex = typename TotalDomain::discrete_element_type;
        using BatchIndex = typename BatchDomain::discrete_element_type;

        static_assert(
                std::is_invocable_v<IntegratorFunction, TotalIndex>,
                "The object passed to Quadrature::operator() is not defined on the total domain.");

        // Get domains
        QuadratureDomain quad_domain(m_coefficients.domain());
        BatchDomain batch_domain(result.domain());

        ddc::ChunkSpan const coeff_proxy = m_coefficients;
        // Loop over batch dimensions
        Kokkos::parallel_for(
                Kokkos::TeamPolicy<>(exec_space, batch_domain.size(), Kokkos::AUTO),
                KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                    const int idx = team.league_rank();
                    BatchIndex ib = to_discrete_element(idx, batch_domain);

                    // Sum over quadrature dimensions
                    double teamSum = 0;
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, quad_domain.size()),
                            [&](int const& thread_index, double& sum) {
                                QuadratureIndex iq = to_discrete_element(thread_index, quad_domain);
                                TotalIndex it(ib, iq);
                                sum += coeff_proxy(iq) * integrated_function(it);
                            },
                            teamSum);
                    result(ib) = teamSum;
                });
    }

private:
    /**
     * A function which converts an integer into a DiscreteElement found in a domain
     * starting from the front. This is useful for iterating over a domain using Kokkos
     * loops.
     *
     * @param[in] idx The index of the requested element.
     * @param[in] dom The domain being iterated over.
     *
     * @return The vector displacement from the front of the domain
     */
    template <class HeadDim, class... DDim>
    KOKKOS_FUNCTION static ddc::DiscreteElement<HeadDim, DDim...> to_discrete_element(
            int idx,
            ddc::DiscreteDomain<HeadDim, DDim...> dom)
    {
        ddc::DiscreteDomain<DDim...> subdomain(dom);
        ddc::DiscreteElement<HeadDim> head_idx(
                ddc::select<HeadDim>(dom).front() + idx / subdomain.size());
        if constexpr (sizeof...(DDim) == 0) {
            return head_idx;
        } else {
            ddc::DiscreteElement<DDim...> tail_idx
                    = to_discrete_element(idx % subdomain.size(), subdomain);
            return ddc::DiscreteElement<HeadDim, DDim...>(head_idx, tail_idx);
        }
    }
};

namespace detail {
template <class NewMemorySpace, class QuadratureDomain, class TotalDomain, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, Quadrature<QuadratureDomain, TotalDomain, MemorySpace>>
{
    using type = Quadrature<QuadratureDomain, TotalDomain, NewMemorySpace>;
};
} // namespace detail
