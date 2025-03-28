

# File rk4.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**rk4.hpp**](rk4_8hpp.md)

[Go to the documentation of this file](rk4_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "vector_field_common.hpp"


template <
        class FieldMem,
        class DerivFieldMem = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class RK4 : public ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>
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

public:
    using base_type::update;

public:
    explicit RK4(IdxRange idx_range) : m_idx_range(idx_range) {}

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
        base_type::copy(y_prime, get_const_field(y));

        // --------- Calculate k1 ------------
        // k1 = f(y)
        dy_calculator(k1, get_const_field(y));

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(y_prime, get_const_field(k1), 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy_calculator(k2, get_const_field(y_prime));

        // --------- Calculate k3 ------------
        // Collect initial conditions
        base_type::copy(y_prime, get_const_field(y));

        // Calculate y_new := y_n + h/2*k_2
        y_update(y_prime, get_const_field(k2), 0.5 * dt);

        // Calculate k3 = f(y_new)
        dy_calculator(k3, get_const_field(y_prime));

        // --------- Calculate k3 ------------
        // Collect initial conditions
        base_type::copy(y_prime, get_const_field(y));

        // Calculate y_new := y_n + h*k_3
        y_update(y_prime, get_const_field(k3), dt);

        // Calculate k4 = f(y_new)
        dy_calculator(k4, get_const_field(y_prime));

        // --------- Update y ------------
        // Calculation of step
        // k_total = k1 + 4 * k2 + k3
        using element_type = typename DerivField::element_type;
        base_type::assemble_k_total(
                exec_space,
                k_total,
                KOKKOS_LAMBDA(std::array<element_type, 4> k) {
                    return k[0] + 2 * k[1] + 2 * k[2] + k[3];
                },
                k1,
                k2,
                k3,
                k4);

        // Calculate y_{n+1} := y_n + (k1 + 2 * k2 + 2 * k3 + k4) * h/6
        y_update(y, get_const_field(k_total), dt / 6.);
    }
};
```


