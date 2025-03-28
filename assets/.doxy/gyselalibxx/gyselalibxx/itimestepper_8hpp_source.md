

# File itimestepper.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**itimestepper.hpp**](itimestepper_8hpp.md)

[Go to the documentation of this file](itimestepper_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include <ddc/ddc.hpp>

#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"

template <
        class FieldMem,
        class DerivFieldMem = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class ITimeStepper
{
    static_assert(
            (ddc::is_chunk_v<FieldMem>) or (is_vector_field_v<FieldMem>)
            or (is_multipatch_field_mem_v<FieldMem>));
    static_assert(
            (ddc::is_chunk_v<DerivFieldMem>) or (is_vector_field_v<DerivFieldMem>)
            or (is_multipatch_field_mem_v<DerivFieldMem>));

    static_assert(
            (std::is_same_v<
                    typename FieldMem::discrete_domain_type,
                    typename DerivFieldMem::discrete_domain_type>)
            || (is_multipatch_field_mem_v<FieldMem> && is_multipatch_field_mem_v<DerivFieldMem>));

    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FieldMem::memory_space>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename DerivFieldMem::memory_space>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");

public:
    using IdxRange = typename FieldMem::discrete_domain_type;


    using ValField = typename FieldMem::span_type;

    using ValConstField = typename FieldMem::view_type;

    using DerivField = typename DerivFieldMem::span_type;

    using DerivConstField = typename DerivFieldMem::view_type;

public:
    void update(ValField y, double dt, std::function<void(DerivField, ValConstField)> dy_calculator)
            const
    {
        update(ExecSpace(), y, dt, dy_calculator);
    }

    void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy_calculator) const
    {
        static_assert(ddc::is_chunk_v<FieldMem>);
        using Idx = typename IdxRange::discrete_element_type;
        update(exec_space, y, dt, dy_calculator, [&](ValField y, DerivConstField dy, double dt) {
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(y),
                    KOKKOS_LAMBDA(Idx const idx) { y(idx) = y(idx) + dy(idx) * dt; });
        });
    }

    virtual void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy_calculator,
            std::function<void(ValField, DerivConstField, double)> y_update) const = 0;

protected:
    void copy(ValField copy_to, ValConstField copy_from) const
    {
        if constexpr (ddc::is_chunk_v<ValField>) {
            ddc::parallel_deepcopy(copy_to, copy_from);
        } else {
            ddcHelper::deepcopy(copy_to, copy_from);
        }
    }

    template <class DerivFieldType, class Idx, class... DDims>
    KOKKOS_FUNCTION static void fill_k_total(
            DerivFieldType k_total,
            Idx i,
            DVector<DDims...> new_val)
    {
        static_assert(
                (std::is_same_v<DerivField, DerivFieldType>)
                || (is_multipatch_field_v<DerivField>));
        ((ddcHelper::get<DDims>(k_total)(i) = ddcHelper::get<DDims>(new_val)), ...);
    }

    template <class FuncType, class... T>
    void assemble_k_total(ExecSpace const& exec_space, DerivField k_total, FuncType func, T... k)
            const
    {
        static_assert(std::conjunction_v<std::is_same<T, DerivField>...>);
        std::size_t constexpr n_args = sizeof...(T);
        using element_type = typename DerivField::element_type;
        static_assert(
                std::is_invocable_r_v<element_type, FuncType, std::array<element_type, n_args>>);
        std::array<DerivField, n_args> k_arr({k...});
        if constexpr (is_vector_field_v<DerivField>) {
            assemble_vector_field_k_total(exec_space, k_total, func, k_arr);
        } else if constexpr (is_multipatch_field_v<DerivField>) {
            assemble_multipatch_field_k_total(exec_space, k_total, func, k_arr);
        } else {
            assemble_field_k_total(exec_space, k_total, func, k_arr);
        }
    }

public:
    template <class FieldType, class FuncType, std::size_t n_args>
    void assemble_field_k_total(
            ExecSpace const& exec_space,
            FieldType k_total,
            FuncType func,
            std::array<FieldType, n_args> k_arr) const
    {
        static_assert(ddc::is_chunk_v<FieldType>);
        using Idx = typename FieldType::discrete_domain_type::discrete_element_type;
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(k_total),
                KOKKOS_LAMBDA(Idx const i) {
                    std::array<double, n_args> k_elems;
                    for (int j(0); j < n_args; ++j) {
                        k_elems[j] = k_arr[j](i);
                    }
                    k_total(i) = func(k_elems);
                });
    }

    template <class FieldType, class FuncType, std::size_t n_args>
    void assemble_vector_field_k_total(
            ExecSpace const& exec_space,
            FieldType k_total,
            FuncType func,
            std::array<FieldType, n_args> k_arr) const
    {
        static_assert(is_vector_field_v<FieldType>);
        using element_type = typename FieldType::element_type;
        using Idx = typename FieldType::discrete_domain_type::discrete_element_type;
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(k_total),
                KOKKOS_LAMBDA(Idx const i) {
                    std::array<element_type, n_args> k_elems;
                    for (int j(0); j < n_args; ++j) {
                        k_elems[j] = k_arr[j](i);
                    }
                    fill_k_total(k_total, i, func(k_elems));
                });
    }

private:
    template <
            class Patch,
            template <typename P>
            typename T,
            class... Patches,
            class FuncType,
            std::size_t n_args>
    void assemble_multipatch_field_k_total_on_patch(
            ExecSpace const& exec_space,
            MultipatchField<T, Patches...> k_total,
            FuncType func,
            std::array<MultipatchField<T, Patches...>, n_args> k_arr) const
    {
        using FieldType = T<Patch>;
        static_assert((ddc::is_chunk_v<FieldType>) or (is_vector_field_v<FieldType>));
        std::array<FieldType, n_args> k_arr_on_patch;
        FieldType k_total_on_patch = k_total.template get<Patch>();
        for (std::size_t i(0); i < n_args; ++i) {
            k_arr_on_patch[i] = k_arr[i].template get<Patch>();
        }
        if constexpr (is_vector_field_v<FieldType>) {
            assemble_vector_field_k_total(exec_space, k_total_on_patch, func, k_arr_on_patch);
        } else {
            assemble_field_k_total(exec_space, k_total_on_patch, func, k_arr_on_patch);
        }
    }

    template <
            template <typename P>
            typename T,
            class... Patches,
            class FuncType,
            std::size_t n_args>
    void assemble_multipatch_field_k_total(
            ExecSpace const& exec_space,
            MultipatchField<T, Patches...> k_total,
            FuncType func,
            std::array<MultipatchField<T, Patches...>, n_args> k_arr) const
    {
        ((assemble_multipatch_field_k_total_on_patch<Patches>(exec_space, k_total, func, k_arr)),
         ...);
    }
};
```


