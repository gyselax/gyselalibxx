

# File crank\_nicolson.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**crank\_nicolson.hpp**](crank__nicolson_8hpp.md)

[Go to the documentation of this file](crank__nicolson_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "l_norm_tools.hpp"
#include "multipatch_math_tools.hpp"
#include "vector_field_common.hpp"


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
    explicit CrankNicolson(
            IdxRange idx_range,
            int const counter = int(20),
            double const epsilon = 1e-12)
        : m_idx_range(idx_range)
        , m_max_counter(counter)
        , m_epsilon(epsilon)
    {
    }

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

            // Re-initialise the characteristic feet
            base_type::copy(y, get_const_field(y_init));

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, get_const_field(k_total), 0.5 * dt);


            // Check convergence
            not_converged
                    = not have_converged(exec_space, get_const_field(y_old), get_const_field(y));


        } while (not_converged and (counter < m_max_counter));
    }


    bool have_converged(ExecSpace const& exec_space, ValConstField y_old, ValConstField y_new) const
    {
        double norm_old = norm_inf(exec_space, y_old);

        double max_diff = error_norm_inf(exec_space, y_old, y_new);

        return (max_diff / norm_old) < m_epsilon;
    }
};
```


