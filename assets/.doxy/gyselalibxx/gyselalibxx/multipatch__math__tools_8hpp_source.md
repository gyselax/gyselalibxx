

# File multipatch\_math\_tools.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**utils**](dir_573def5310cd01d120c251a7885d602c.md) **>** [**multipatch\_math\_tools.hpp**](multipatch__math__tools_8hpp.md)

[Go to the documentation of this file](multipatch__math__tools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "l_norm_tools.hpp"
#include "multipatch_field.hpp"

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
```


