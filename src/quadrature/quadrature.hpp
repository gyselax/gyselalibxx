// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"

/**
 * @brief A class providing an operator for integrating functions defined on a discrete domain.
 */
template <class QuadratureDomain, class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space>
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
};

namespace detail {
template <class NewMemorySpace, class QuadratureDomain, class MemorySpace>
struct OnMemorySpace<NewMemorySpace, Quadrature<QuadratureDomain, MemorySpace>>
{
    using type = Quadrature<QuadratureDomain, NewMemorySpace>;
};
} // namespace detail
