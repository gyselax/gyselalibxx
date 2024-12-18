// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "math_tools.hpp"
#include "multipatch_math_tools.hpp"
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
        using element_type = typename DerivField::element_type;

        FieldMem y_init_alloc(m_idx_range);
        FieldMem y_old_alloc(m_idx_range);
        DerivFieldMem k1_alloc(m_idx_range);
        DerivFieldMem k_new_alloc(m_idx_range);
        DerivFieldMem k_total_alloc(m_idx_range);
        ValField y_init = get_field(y_init_alloc);
        ValField y_old = get_field(y_old_alloc);
        DerivField k1 = get_field(k1_alloc);
        DerivField k_new = get_field(k_new_alloc);
        DerivField k_total = get_field(k_total_alloc);


        base_type::copy(y_init, get_const_field(y));

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(k1, get_const_field(y));

        // -------- Calculate k_new ----------
        bool not_converged = true;
        int counter = 0;
        do {
            counter++;

            // Calculate k_new = f(y_new)
            dy_calculator(k_new, get_const_field(y));

            // Calculation of step
            // k_total = k1 + k_new
            base_type::assemble_k_total(
                    exec_space,
                    k_total,
                    KOKKOS_LAMBDA(std::array<element_type, 2> k) { return k[0] + k[1]; },
                    k1,
                    k_new);

            // Save the old characteristic feet
            base_type::copy(y_old, get_const_field(y));

            // Re-initialiase the characteristic feet
            base_type::copy(y, get_const_field(y_init));

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, get_const_field(k_total), 0.5 * dt);


            // Check convergence
            not_converged
                    = not have_converged(exec_space, get_const_field(y_old), get_const_field(y));


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
        double norm_old = norm_inf(exec_space, y_old);

        double max_diff = error_norm_inf(exec_space, y_old, y_new);

        return (max_diff / norm_old) < m_epsilon;
    }
};
