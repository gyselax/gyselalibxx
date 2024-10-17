// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "utils_tools.hpp"
#include "vector_field_common.hpp"


/**
 * @brief A class which provides an implementation of a Crank-Nicolson method.
 *
 * A class which provides an implementation of a Crank-Nicolson method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * The values which evolve are defined on an index range.
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the Crank-Nicolson method is given by :
 * @f$ y^{k} =  y^{n} + \frac{dt}{2} \left(f(t^{n}, y^{n}) + f(t^{k}, y^{k}) \right)@f$.
 *
 * The method is an implicit method.
 * If @f$ |y^{k+1} -  y^{k}| < \varepsilon @f$, then we set @f$ y^{n+1} = y^{k+1} @f$.
 *
 * The method is order 2.
 *
 */
template <
        class FieldMem,
        class DerivFieldMem = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class CrankNicolson : public ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>
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
    int const m_max_counter;
    double const m_epsilon;

public:
    using base_type::update;

public:
    /**
     * @brief Create a CrankNicolson object.
     * @param[in] idx_range
     *      The index range on which the points which evolve over time are defined.
     * @param[in] counter
     *      The maximal number of loops for the implicit method.
     * @param[in] epsilon
     *      The @f$ \varepsilon @f$ upperbound of the difference of two steps
     *      in the implicit method: @f$ |y^{k+1} -  y^{k}| < \varepsilon @f$.
     */
    explicit CrankNicolson(
            IdxRange idx_range,
            int const counter = int(20),
            double const epsilon = 1e-12)
        : m_idx_range(idx_range)
        , m_max_counter(counter)
        , m_epsilon(epsilon)
    {
    }

    /**
     * @brief Carry out one step of the Crank-Nicolson scheme.
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
        FieldMem m_y_init_alloc(m_idx_range);
        FieldMem m_y_old_alloc(m_idx_range);
        DerivFieldMem m_k1_alloc(m_idx_range);
        DerivFieldMem m_k_new_alloc(m_idx_range);
        DerivFieldMem m_k_total_alloc(m_idx_range);
        ValField m_y_init = get_field(m_y_init_alloc);
        ValField m_y_old = get_field(m_y_old_alloc);
        DerivField m_k1 = get_field(m_k1_alloc);
        DerivField m_k_new = get_field(m_k_new_alloc);
        DerivField m_k_total = get_field(m_k_total_alloc);


        copy(m_y_init, y);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(m_k1, y);

        // -------- Calculate k_new ----------
        bool not_converged = true;
        int counter = 0;
        do {
            counter++;

            // Calculate k_new = f(y_new)
            dy_calculator(m_k_new, y);

            // Calculation of step
            if constexpr (is_field_v<DerivFieldMem>) {
                ddc::parallel_for_each(
                        exec_space,
                        get_idx_range(m_k_total),
                        KOKKOS_CLASS_LAMBDA(Idx const i) {
                            // k_total = k1 + k_new
                            fill_k_total(i, m_k_total, m_k1(i) + m_k_new(i));
                        });
            } else {
                ddc::parallel_for_each(
                        exec_space,
                        get_idx_range(m_k_total),
                        KOKKOS_CLASS_LAMBDA(Idx const i) {
                            // k_total = k1 + k_new
                            m_k_total(i) = m_k1(i) + m_k_new(i);
                        });
            }

            // Save the old characteristic feet
            copy(m_y_old, y);

            // Re-initiliase the characteristic feet
            copy(y, m_y_init);

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, m_k_total, 0.5 * dt);


            // Check convergence
            not_converged = not have_converged(exec_space, m_y_old, y);


        } while (not_converged and (counter < m_max_counter));
    }


    /**
     * Check if the relative difference of the function between
     * two time steps is below epsilon.
     *
     * This function should be private. It is not due to the inclusion of a KOKKOS_LAMBDA
     * function.
     *
     * @param[in] exec_space
     *     The space on which the function is executed (CPU/GPU).
     * @param[in] y_old
     *     The value of the function at the previous time step.
     * @param[in] y_new
     *     The updated value of the function at the new time step.
     * @returns
     *     True if converged, False otherwise.
     */
    bool have_converged(ExecSpace const& exec_space, ValConstField y_old, ValConstField y_new) const
    {
        IdxRange const idx_range = get_idx_range(y_old);

        double norm_old = ddc::parallel_transform_reduce(
                exec_space,
                idx_range,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(Idx const idx) { return norm_inf(y_old(idx)); });

        double max_diff = ddc::parallel_transform_reduce(
                exec_space,
                idx_range,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(Idx const idx) { return norm_inf(y_old(idx) - y_new(idx)); });

        return (max_diff / norm_old) < m_epsilon;
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
    KOKKOS_FUNCTION void fill_k_total(Idx i, DerivField m_k_total, Coord<DDims...> new_val) const
    {
        ((ddcHelper::get<DDims>(m_k_total)(i) = ddc::get<DDims>(new_val)), ...);
    }
};
