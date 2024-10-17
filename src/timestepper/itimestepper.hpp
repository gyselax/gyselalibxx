// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include <ddc/ddc.hpp>

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
    static_assert(ddc::is_chunk_v<FieldMem> or is_field_v<FieldMem>);
    static_assert(ddc::is_chunk_v<DerivFieldMem> or is_field_v<DerivFieldMem>);

    static_assert(std::is_same_v<
                  typename FieldMem::discrete_domain_type,
                  typename DerivFieldMem::discrete_domain_type>);

    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename FieldMem::memory_space>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, typename DerivFieldMem::memory_space>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");

public:
    /// The type of the index range on which the values of the function are defined.
    using IdxRange = typename FieldMem::discrete_domain_type;

    /// The type of the indices of the function.
    using Idx = typename IdxRange::discrete_element_type;

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
};
