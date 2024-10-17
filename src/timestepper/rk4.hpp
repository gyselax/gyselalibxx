// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "vector_field_common.hpp"


/**
 * @brief A class which provides an implementation of a fourth-order Runge-Kutta method.
 *
 * A class which provides an implementation of a fourth-order Runge-Kutta method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * The values which evolve are defined on an index range.
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the Runge-Kutta 3 method is given by :
 * @f$ y^{n+1} =  y^{n} + \frac{dt}{6} \left(k_1 + 2k_2 + 2k_3 + k_4 \right) @f$,
 *
 * with
 *
 * - @f$ k_1 = f(t^{n}, y^{n}) @f$,
 * - @f$ k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 ) @f$,
 * - @f$ k_3 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_2 ) @f$,
 * - @f$ k_3 = f(t^{n}, y^{n} + dt k_3 ) @f$.
 */
template <
        class FieldMem,
        class DerivFieldMem = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class RK4 : public ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>
{
    using base_type = ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>;

public:
    using typename base_type::IdxRange;

    using typename base_type::Idx;

    using typename base_type::ValConstField;
    using typename base_type::ValField;

    using typename base_type::DerivConstField;
    using typename base_type::DerivField;

private:
    IdxRange const m_idx_range;

public:
    using base_type::update;

public:
    /**
     * @brief Create a RK4 object.
     * @param[in] idx_range The index range on which the points which evolve over time are defined.
     */
    explicit RK4(IdxRange idx_range) : m_idx_range(idx_range) {}

    /**
     * @brief Carry out one step of the Runge-Kutta scheme.
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
    void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy_calculator,
            std::function<void(ValField, DerivConstField, double)> y_update) const final
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename FieldMem::memory_space>::accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename DerivFieldMem::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        FieldMem y_prime_alloc(m_idx_range);
        DerivFieldMem k1_alloc(m_idx_range);
        DerivFieldMem k2_alloc(m_idx_range);
        DerivFieldMem k3_alloc(m_idx_range);
        DerivFieldMem k4_alloc(m_idx_range);
        DerivFieldMem k_total_alloc(m_idx_range);

        ValField y_prime = get_field(y_prime_alloc);
        DerivField k1 = get_field(k1_alloc);
        DerivField k2 = get_field(k2_alloc);
        DerivField k3 = get_field(k3_alloc);
        DerivField k4 = get_field(k4_alloc);
        DerivField k_total = get_field(k_total_alloc);


        // Save initial conditions
        copy(y_prime, y);

        // --------- Calculate k1 ------------
        // k1 = f(y)
        dy_calculator(k1, y);

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(y_prime, k1, 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy_calculator(k2, y_prime);

        // --------- Calculate k3 ------------
        // Collect initial conditions
        copy(y_prime, y);

        // Calculate y_new := y_n + h/2*k_2
        y_update(y_prime, k2, 0.5 * dt);

        // Calculate k3 = f(y_new)
        dy_calculator(k3, y_prime);

        // --------- Calculate k3 ------------
        // Collect initial conditions
        copy(y_prime, y);

        // Calculate y_new := y_n + h*k_3
        y_update(y_prime, k3, dt);

        // Calculate k4 = f(y_new)
        dy_calculator(k4, y_prime);

        // --------- Update y ------------
        // Calculation of step
        if constexpr (is_field_v<DerivFieldMem>) {
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(k_total),
                    KOKKOS_CLASS_LAMBDA(Idx const i) {
                        // k_total = k1 + 4 * k2 + k3
                        fill_k_total(i, k_total, k1(i) + 2 * k2(i) + 2 * k3(i) + k4(i));
                    });
        } else {
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(k_total),
                    KOKKOS_LAMBDA(Idx const i) {
                        // k_total = k1 + 4 * k2 + k3
                        k_total(i) = k1(i) + 2 * k2(i) + 2 * k3(i) + k4(i);
                    });
        }

        // Calculate y_{n+1} := y_n + (k1 + 2 * k2 + 2 * k3 + k4) * h/6
        y_update(y, k_total, dt / 6.);
    }

private:
    void copy(ValField copy_to, ValConstField copy_from) const
    {
        if constexpr (is_field_v<ValField>) {
            ddcHelper::deepcopy(copy_to, copy_from);
        } else {
            ddc::parallel_deepcopy(copy_to, copy_from);
        }
    }

    template <class... DDims>
    KOKKOS_FUNCTION void fill_k_total(Idx i, DerivField k_total, Coord<DDims...> new_val) const
    {
        ((ddcHelper::get<DDims>(k_total)(i) = ddc::get<DDims>(new_val)), ...);
    }
};
