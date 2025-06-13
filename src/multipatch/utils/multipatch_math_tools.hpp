// SPDX-License-Identifier: MIT
#pragma once
#include "l_norm_tools.hpp"
#include "multipatch_field.hpp"

/**
 * @brief Compute the infinity norm for a Field or VectorField over multiple patches.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] multipatch_function The function whose norm is calculated.
 * @return A double containing the value of the infinity norm.
 */
template <class ExecSpace, template <typename P> typename T, class... Patches>
double norm_inf(ExecSpace exec_space, MultipatchField<T, Patches...> multipatch_function)
{
    using FuncType = MultipatchField<T, Patches...>;
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    constexpr std::size_t NPatches = multipatch_function.size();
    std::array<double, NPatches> norm_inf_on_patch(
            {(norm_inf(exec_space, multipatch_function.template get<Patches>()))...});
    double result(0.0);
    for (std::size_t i(0); i < NPatches; ++i) {
        result = std::max(result, norm_inf_on_patch[i]);
    }
    return result;
}

/**
 * @brief Compute the infinity norm of the error between 2 Fields or VectorFields over multiple patches.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] multipatch_function The calculated function.
 * @param[in] multipatch_exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinity norm.
 */
template <class ExecSpace, template <typename P> typename T, class... Patches>
double error_norm_inf(
        ExecSpace exec_space,
        MultipatchField<T, Patches...> multipatch_function,
        MultipatchField<T, Patches...> multipatch_exact_function)
{
    using FuncType = MultipatchField<T, Patches...>;
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    constexpr std::size_t NPatches = multipatch_function.size();
    std::array<double, NPatches> norm_inf_on_patch({(error_norm_inf(
            exec_space,
            multipatch_function.template get<Patches>(),
            multipatch_exact_function.template get<Patches>()))...});
    double result(0.0);
    for (std::size_t i(0); i < NPatches; ++i) {
        result = std::max(result, norm_inf_on_patch[i]);
    }
    return result;
}
