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
template <class ExecSpace, class... Patches, class TypeSeqInternalTypes>
double norm_inf(
        ExecSpace exec_space,
        detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes>
                multipatch_function)
{
    using FuncType
            = detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes>;
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
template <class ExecSpace, class... Patches, class TypeSeqInternalTypes>
double error_norm_inf(
        ExecSpace exec_space,
        detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes>
                multipatch_function,
        detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes>
                multipatch_exact_function)
{
    using FuncType
            = detail::MultipatchField<ddc::detail::TypeSeq<Patches...>, TypeSeqInternalTypes>;
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
