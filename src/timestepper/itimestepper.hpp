// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include <ddc/ddc.hpp>

#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"

/**
 * @brief The superclass from which all timestepping methods inherit.
 *
 * The class exposes three update functions which are used to carry out one step
 * of the chosen timestepping method to solve an ODE of the form:
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 */
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
    /// The type of the index range on which the values of the function are defined.
    using IdxRange = typename FieldMem::discrete_domain_type;


    /// The type of the values of the function being evolved.
    using ValField = typename FieldMem::span_type;

    /// The constant type of the values of the function being evolved.
    using ValConstField = typename FieldMem::view_type;

    /// The type of the derivatives of the function being evolved.
    using DerivField = typename DerivFieldMem::span_type;

    /// The constant type of the derivatives values of the function being evolved.
    using DerivConstField = typename DerivFieldMem::view_type;

public:
    /**
     * @brief Carry out one step of the timestepping scheme.
     *
     * This function is a wrapper around the update function below. The values of the function are
     * updated using the trivial method $f += df * dt$. This is the standard method however some
     * cases may need a more complex update function which is why the more explicit method is
     * also provided.
     *
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the index range.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy_calculator
     *     The function describing how the derivative of the evolve function is calculated.
     */
    void update(ValField y, double dt, std::function<void(DerivField, ValConstField)> dy_calculator)
            const
    {
        update(ExecSpace(), y, dt, dy_calculator);
    }

    /**
     * @brief Carry out one step of the timestepping scheme.
     *
     * This function is a wrapper around the update function below. The values of the function are
     * updated using the trivial method $f += df * dt$. This is the standard method however some
     * cases may need a more complex update function which is why the more explicit method is
     * also provided.
     *
     * @param[in] exec_space
     *     The space on which the function is executed (CPU/GPU).
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the index range.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy_calculator
     *     The function describing how the derivative of the evolve function is calculated.
     */
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

    /**
     * @brief Carry out one step of the timestepping scheme.
     *
     * @param[in] exec_space
     *     The space on which the function is executed (CPU/GPU).
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the index range.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy_calculator
     *     The function describing how the derivative of the evolve function is calculated.
     * @param[in] y_update
     *     The function describing how the value(s) are updated using the derivative.
     */
    virtual void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy_calculator,
            std::function<void(ValField, DerivConstField, double)> y_update) const = 0;

protected:
    /**
     * @brief Make a copy of the values of the function being evolved.
     *
     * @param[out] copy_to the field that the values should be copied to.
     * @param[in] copy_from The field that the values should be copied from.
     */
    void copy(ValField copy_to, ValConstField copy_from) const
    {
        if constexpr (ddc::is_chunk_v<ValField>) {
            ddc::parallel_deepcopy(copy_to, copy_from);
        } else {
            ddcHelper::deepcopy(copy_to, copy_from);
        }
    }

    /**
     * @brief A method to fill an element of a vector field.
     *
     * @param[out] k_total The vector field that will be filled.
     * @param[in] i The index where the vector field should be filled.
     * @param[in] new_val The coordinate that should be saved to the vector field.
     */
    template <class DerivFieldType, class Idx, class... DDims>
    KOKKOS_FUNCTION static void fill_k_total(DerivFieldType k_total, Idx i, Coord<DDims...> new_val)
    {
        static_assert(
                (std::is_same_v<DerivField, DerivFieldType>)
                || (is_multipatch_field_v<DerivField>));
        ((ddcHelper::get<DDims>(k_total)(i) = ddc::get<DDims>(new_val)), ...);
    }

    /**
     * @brief A method to assemble multiple derivative fields into one. This method is responsible
     * for choosing how this is done depending on the type of the derivative field.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k The derivative fields being combined.
     */
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
    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a Field (ddc::ChunkSpan).
     * This function should be private but is public due to Cuda restrictions.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
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

    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a VectorField.
     * This function should be private but is public due to Cuda restrictions.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
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
    /**
     * Calculate func(k_arr[0], k_arr[1], ...) on one patch of a MultipatchField.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
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

    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a MultipatchField.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
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
