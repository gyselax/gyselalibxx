

# File rk2.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**rk2.hpp**](rk2_8hpp.md)

[Go to the documentation of this file](rk2_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "vector_field_common.hpp"

template <class ValType, class DerivType = ValType, class ExecSpace = Kokkos::DefaultExecutionSpace>
class RK2
{
    static_assert(!timestepper_detail::FieldLike<ValType>);
    static_assert(!timestepper_detail::FieldLike<DerivType>);

public:
    using ValFieldMem = ValType;

    using DerivFieldMem = DerivType;

    using exec_space = ExecSpace;

public:
    explicit KOKKOS_DEFAULTED_FUNCTION RK2() = default;

    template <
            class DYFunctor,
            class YFunctor
            = decltype(timestepper_detail::serial_y_update<ValType&, DerivType const&>)>
    KOKKOS_FUNCTION void update(
            ValType& y,
            double dt,
            DYFunctor dy_calculator,
            YFunctor y_update
            = timestepper_detail::serial_y_update<ValType&, DerivType const&>) const
    {
        static_assert(std::is_invocable_v<DYFunctor, DerivType&, ValType>);
        ValType y_prime;
        DerivType k1;
        DerivType k2;

        // Save initial conditions
        timestepper_detail::copy_helper<ValType>::copy(y_prime, y);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y)
        dy_calculator(k1, y);

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(y_prime, k1, 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy_calculator(k2, y_prime);

        // ----------- Update y --------------
        // Calculate y_{n+1} := y_n + h*k_2
        y_update(y, k2, dt);
    }
};

template <
        timestepper_detail::FieldLike FieldMem,
        timestepper_detail::FieldLike DerivFieldMem,
        class ExecSpace>
class RK2<FieldMem, DerivFieldMem, ExecSpace>
    : public ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>
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
    explicit RK2(IdxRange idx_range) : m_idx_range(idx_range) {}

    void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy_calculator,
            std::function<void(ValField, DerivConstField, double)> y_update) const final
    {
        DerivFieldMem k1_alloc("k1 (RK2::update)", m_idx_range);
        DerivFieldMem k2_alloc("k2 (RK2::update)", m_idx_range);
        FieldMem y_prime_alloc("y_prime (RK2::update)", m_idx_range);

        DerivField k1 = get_field(k1_alloc);
        DerivField k2 = get_field(k2_alloc);
        ValField y_prime = get_field(y_prime_alloc);

        // Save initial conditions
        timestepper_detail::copy_helper<FieldMem>::copy(y_prime, get_const_field(y));

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y)
        dy_calculator(k1, get_const_field(y));

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(y_prime, get_const_field(k1), 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy_calculator(k2, get_const_field(y_prime));

        // ----------- Update y --------------
        // Calculate y_{n+1} := y_n + h*k_2
        y_update(y, get_const_field(k2), dt);
    }
};

using RK2Builder = ExplicitTimeStepperBuilder<RK2>;
```


