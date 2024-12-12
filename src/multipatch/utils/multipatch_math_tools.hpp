// SPDX-License-Identifier: MIT
#pragma once
#include "Lnorm_tools.hpp"
#include "multipatch_field.hpp"

/**
 * @brief Compute the infinity norm for a Field or VectorField over multiple patches.
 * @param[in] exec_space The space on which the function is executed (CPU/GPU).
 * @param[in] function The function whose norm is calcuated.
 * @return A double containing the value of the infinty norm.
 */
template <class ExecSpace, template <typename P> typename T, class... Patches>
double norm_inf(ExecSpace exec_space, MultipatchField<T, Patches...> multipatch_function)
{
    using FuncType = MultipatchField<T, Patches...>;
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FuncType::memory_space>::accessible);
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    constexpr std::size_t NPatches = multipatch_function.size();
    IdxRangeFunc idx_range = get_idx_range(multipatch_function);
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
 * @param[in] function The calculated function.
 * @param[in] exact_function The exact function with which the calculated function is compared.
 * @return A double containing the value of the infinty norm.
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
    using IdxRangeFunc = typename FuncType::discrete_domain_type;
    constexpr std::size_t NPatches = multipatch_function.size();
    IdxRangeFunc idx_range = get_idx_range(multipatch_function);
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
