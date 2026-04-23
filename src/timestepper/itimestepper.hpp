// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"
#include "tensor.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

/// @cond
namespace timestepper_detail {

template <typename T>
concept FieldLike = requires
{
    typename T::span_type;
    typename T::view_type;
};

template <class T>
struct GetReferenceField
{
    static_assert(!FieldLike<T>);
    using type = T&;
};

template <class T>
struct GetConstReferenceField
{
    static_assert(!FieldLike<T>);
    using type = T const&;
};

template <FieldLike T>
struct GetReferenceField<T>
{
    using type = typename T::span_type;
};

template <FieldLike T>
struct GetConstReferenceField<T>
{
    using type = typename T::view_type;
};

template <class T>
using reference_t = typename GetReferenceField<T>::type;
template <class T>
using const_reference_t = typename GetConstReferenceField<T>::type;

template <class FieldMem, class DerivFieldMemType>
concept Compatible
        = !(FieldLike<FieldMem> && FieldLike<DerivFieldMemType>)
          || (std::is_same_v<
                  typename FieldMem::discrete_domain_type,
                  typename DerivFieldMemType::discrete_domain_type>)
          || (is_multipatch_field_mem_v<FieldMem> && is_multipatch_field_mem_v<DerivFieldMemType>);

template <class ExecSpace, class FieldMem>
concept Accessible
        = !(FieldLike<FieldMem>)
          || bool(Kokkos::SpaceAccessibility<ExecSpace, typename FieldMem::memory_space>::
                          accessible);

template <class FieldMem>
concept ExpectedTimeStepperType
        = (ddc::is_chunk_v<FieldMem>) || (is_vector_field_v<FieldMem>)
          || (is_multipatch_field_mem_v<FieldMem>) || (std::is_floating_point_v<FieldMem>)
          || (is_tensor_type_v<FieldMem>) || (ddcHelper::is_coordinate_v<FieldMem>);

template <class T>
struct IdxRangeType
{
    static_assert(!FieldLike<T>);
    using type = IdxRange<>;
};

template <FieldLike T>
struct IdxRangeType<T>
{
    using type = typename T::discrete_domain_type;
};

template <class T>
struct ElementType
{
    static_assert(!FieldLike<T>);
    using type = std::remove_reference_t<T>;
};

template <FieldLike T>
struct ElementType<T>
{
    using type = typename T::element_type;
};

template <class FieldType>
struct copy_helper
{
    static_assert(!FieldLike<FieldType>);
    /**
     * @brief Make a copy of the values of the function being evolved.
     *
     * @param[out] copy_to the field that the values should be copied to.
     * @param[in] copy_from The field that the values should be copied from.
     */
    static KOKKOS_FUNCTION void copy(FieldType& copy_to, FieldType const& copy_from)
    {
        copy_to = copy_from;
    }
};

template <timestepper_detail::FieldLike FieldMemType>
struct copy_helper<FieldMemType>
{
    static_assert(!ddc::is_chunk_v<FieldMemType>);
    static void copy(
            reference_t<FieldMemType> copy_to,
            const_reference_t<FieldMemType> const& copy_from)
    {
        ddcHelper::deepcopy(copy_to, copy_from);
    }
};

template <class ElementType, class IdxRangeType, class MemSpace>
struct copy_helper<FieldMem<ElementType, IdxRangeType, MemSpace>>
{
    static void copy(
            Field<ElementType, IdxRangeType, MemSpace> copy_to,
            ConstField<ElementType, IdxRangeType, MemSpace> copy_from)
    {
        ddc::parallel_deepcopy(copy_to, copy_from);
    }
};

/**
 * @brief A method to fill an element of a vector field.
 *
 * @param[out] k_total The vector field that will be filled.
 * @param[in] i The index where the vector field should be filled.
 * @param[in] new_val The coordinate that should be saved to the vector field.
 */
template <class DerivFieldType, class Idx, class... DDims>
KOKKOS_FUNCTION void fill_k_total(DerivFieldType k_total, Idx i, DVector<DDims...> new_val)
{
    ((ddcHelper::get<DDims>(k_total)(i) = ddcHelper::get<DDims>(new_val)), ...);
}

template <class ExecSpace, class DerivFieldType>
struct assemble_helper
{
    static_assert(!FieldLike<DerivFieldType>);

    template <class FuncType, class... T>
    static KOKKOS_FUNCTION void assemble_k_total(DerivFieldType& k_total, FuncType func, T... k)
    {
        std::array<DerivFieldType, sizeof...(T)> k_arr({k...});
        k_total = func(k_arr);
    }
};

template <class ExecSpace, class ElementType, class IdxRangeType, class MemSpace>
struct assemble_helper<ExecSpace, FieldMem<ElementType, IdxRangeType, MemSpace>>
{
    using DerivFieldType = Field<ElementType, IdxRangeType, MemSpace>;

    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a Field (ddc::ChunkSpan).
     * This function should be private but is public due to Cuda restrictions.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
    template <class FuncType, class... T>
    static void assemble_k_total(
            ExecSpace const& exec_space,
            DerivFieldType k_total,
            FuncType func,
            T... k)
    {
        std::size_t constexpr n_args = sizeof...(T);
        std::array<DerivFieldType, n_args> k_arr({k...});
        using Idx = typename IdxRangeType::discrete_element_type;
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                exec_space,
                get_idx_range(k_total),
                KOKKOS_LAMBDA(Idx const i) {
                    std::array<ElementType, n_args> k_elems;
                    for (int j(0); j < n_args; ++j) {
                        k_elems[j] = k_arr[j](i);
                    }
                    k_total(i) = func(k_elems);
                });
    }
};

template <
        class ExecSpace,
        class ElementType,
        class IdxRangeType,
        class VectorIndexSet,
        class MemSpace>
struct assemble_helper<
        ExecSpace,
        VectorFieldMem<ElementType, IdxRangeType, VectorIndexSet, MemSpace>>
{
    using DerivFieldType = VectorField<ElementType, IdxRangeType, VectorIndexSet, MemSpace>;

    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a VectorField.
     * This function should be private but is public due to Cuda restrictions.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
    template <class FuncType, class... T>
    static void assemble_k_total(
            ExecSpace const& exec_space,
            DerivFieldType k_total,
            FuncType func,
            T... k)
    {
        std::size_t constexpr n_args = sizeof...(T);
        std::array<DerivFieldType, n_args> k_arr({k...});
        using element_type = Tensor<ElementType, VectorIndexSet>;
        using Idx = typename IdxRangeType::discrete_element_type;
        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
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
};

template <class ExecSpace, template <typename P> typename T, class... Patches>
struct assemble_helper<ExecSpace, MultipatchFieldMem<T, Patches...>>
{
    using DerivFieldType = MultipatchFieldMem<T, Patches...>::span_type;

    /**
     * Calculate func(k_arr[0], k_arr[1], ...) when FieldType is a MultipatchField.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] k_total The field to be filled with the combined derivative fields.
     * @param[in] func A function which combines an element from each of the derivative fields.
     * @param[in] k_arr The derivative fields being combined.
     */
    template <class FuncType, class... KType>
    static void assemble_k_total(
            ExecSpace const& exec_space,
            DerivFieldType k_total,
            FuncType func,
            KType... k)
    {
        ((assemble_multipatch_field_k_total_on_patch<Patches>(exec_space, k_total, func, k...)),
         ...);
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
    template <class Patch, class FuncType, class... KType>
    static void assemble_multipatch_field_k_total_on_patch(
            ExecSpace const& exec_space,
            DerivFieldType k_total,
            FuncType func,
            KType... k)
    {
        timestepper_detail::assemble_helper<ExecSpace, T<Patch>>::assemble_k_total(
                exec_space,
                k_total.template get<Patch>(),
                func,
                k.template get<Patch>()...);
    }
};

template <
        class ValField,
        class DerivConstField,
        typename = std::enable_if_t<!FieldLike<ValField>>,
        typename = std::enable_if_t<!FieldLike<DerivConstField>>>
KOKKOS_FUNCTION void serial_y_update(ValField y, DerivConstField dy, double dt)
{
    y += dy * dt;
}

} // namespace timestepper_detail
/// @endcond

/**
 * @brief The superclass from which all timestepping methods inherit.
 *
 * The class exposes three update functions which are used to carry out one step
 * of the chosen timestepping method to solve an ODE of the form:
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 */
template <
        class FieldMem,
        class DerivFieldMemType = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class ITimeStepper
{
    static_assert(timestepper_detail::FieldLike<FieldMem>);
    static_assert(timestepper_detail::FieldLike<DerivFieldMemType>);
    static_assert(timestepper_detail::Compatible<FieldMem, DerivFieldMemType>);

    static_assert(
            timestepper_detail::Accessible<ExecSpace, FieldMem>,
            "MemorySpace has to be accessible for ExecutionSpace.");
    static_assert(
            timestepper_detail::Accessible<ExecSpace, DerivFieldMemType>,
            "MemorySpace has to be accessible for ExecutionSpace.");

public:
    /// The type of the index range on which the values of the function are defined.
    using IdxRange = typename timestepper_detail::IdxRangeType<FieldMem>::type;

    /// The type of the memory allocation for the values of the function being evolved.
    using ValFieldMem = FieldMem;

    /// The type of the values of the function being evolved.
    using ValField = timestepper_detail::reference_t<FieldMem>;

    /// The constant type of the values of the function being evolved.
    using ValConstField = timestepper_detail::const_reference_t<FieldMem>;

    /// The type of the memory allocation for the derivatives of the function being evolved.
    using DerivFieldMem = DerivFieldMemType;

    /// The type of the derivatives of the function being evolved.
    using DerivField = timestepper_detail::reference_t<DerivFieldMem>;

    /// The constant type of the derivatives values of the function being evolved.
    using DerivConstField = timestepper_detail::const_reference_t<DerivFieldMem>;

    /// The space (CPU/GPU) where the calculations are carried out.
    using exec_space = ExecSpace;

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
        using Idx = typename IdxRange::discrete_element_type;
        update(exec_space, y, dt, dy_calculator, [&](ValField y, DerivConstField dy, double dt) {
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
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
};

/**
 * @brief A class to indicate that an explicit time stepper should be constructed for use in other operators.
 *
 * This class is a time stepper builder. A time stepper builder is designed to construct a
 * time stepper upon request. This allows the simulation to choose the method without
 * needing to know the specifics of the types with which it should be initialised.
 * This class should be specialised for the explicit time stepper builders.
 */
template <template <class FieldMem, class DerivFieldMem, class ExecSpace> typename TimeStepper>
class ExplicitTimeStepperBuilder
{
public:
    /**
     * @brief A constructor for the TimeStepperBuilder
     */
    ExplicitTimeStepperBuilder() {}

    /**
     * The type of the TimeStepper that will be constructed to solve an equation whose field
     * and derivative(s) have the specified type.
     * @tparam FieldMem The type of the data storage for the function.
     * @tparam DerivFieldMem The type of the data storage for the derivative of the function.
     * @tparam ExecSpace The space (CPU/GPU) where the calculations are carried out.
     *                   This template parameter is ignored if the FieldMem is a scalar.
     */
    template <
            class FieldMem,
            class DerivFieldMem = FieldMem,
            class ExecSpace = Kokkos::DefaultExecutionSpace>
    using time_stepper_t = TimeStepper<FieldMem, DerivFieldMem, ExecSpace>;

    /**
     * @brief Allocate the TimeStepper object for FieldLike types.
     * @tparam ChosenTimeStepper The type of the TimeStepper to be constructed (obtained from time_stepper_t).
     * @param[in] idx_range The index range on which the operator will act (and allocate memory).
     */
    template <class ChosenTimeStepper>
    static auto preallocate(typename ChosenTimeStepper::IdxRange const idx_range)
    {
        static_assert(
                timestepper_detail::FieldLike<typename ChosenTimeStepper::ValFieldMem>,
                "An index range should not be provided to preallocate for scalar timesteppers.");
        static_assert(std::is_same_v<
                      ChosenTimeStepper,
                      time_stepper_t<
                              typename ChosenTimeStepper::ValFieldMem,
                              typename ChosenTimeStepper::DerivFieldMem,
                              typename ChosenTimeStepper::exec_space>>);
        return ChosenTimeStepper(idx_range);
    }

    /**
     * @brief Allocate the TimeStepper object for scalar (non-FieldLike) types.
     * @tparam ChosenTimeStepper The type of the TimeStepper to be constructed (obtained from time_stepper_t).
     */
    template <class ChosenTimeStepper>
    static auto preallocate()
    {
        static_assert(
                !timestepper_detail::FieldLike<typename ChosenTimeStepper::ValFieldMem>,
                "An index range must be provided to preallocate for FieldLike timesteppers.");
        static_assert(std::is_same_v<
                      ChosenTimeStepper,
                      time_stepper_t<
                              typename ChosenTimeStepper::ValFieldMem,
                              typename ChosenTimeStepper::DerivFieldMem,
                              typename ChosenTimeStepper::exec_space>>);
        return ChosenTimeStepper();
    }
};

namespace detail {

template <class T>
inline constexpr bool enable_is_timestepper_builder = false;

template <template <class FieldMem, class DerivFieldMem, class ExecSpace> typename TimeStepper>
inline constexpr bool enable_is_timestepper_builder<ExplicitTimeStepperBuilder<TimeStepper>> = true;

} // namespace detail

template <typename Type>
inline constexpr bool is_timestepper_builder_v
        = detail::enable_is_timestepper_builder<std::remove_const_t<std::remove_reference_t<Type>>>;
