

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


template <class ValType, class DerivType = ValType, class ExecSpace = Kokkos::DefaultExecutionSpace>
class CrankNicolson
{
    static_assert(!timestepper_detail::FieldLike<ValType>);
    static_assert(!timestepper_detail::FieldLike<DerivType>);

public:
    using ValFieldMem = ValType;

    using DerivFieldMem = DerivType;

    using exec_space = ExecSpace;

private:
    int const m_max_counter;
    double const m_epsilon;

public:
    explicit CrankNicolson(int const counter = int(20), double const epsilon = 1e-12)
        : m_max_counter(counter)
        , m_epsilon(epsilon)
    {
    }

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
        ValType y_init;
        ValType y_old;
        DerivType k1;
        DerivType k_new;
        DerivType k_total;

        // Save initial conditions
        timestepper_detail::copy_helper<ValType>::copy(y_init, y);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(k1, y);

        // -------- Calculate k_new ----------
        bool not_converged = true;
        int counter = 0;
        do {
            counter++;

            // Calculate k_new = f(y_new)
            dy_calculator(k_new, y);

            // Calculation of step
            // k_total = k1 + k_new
            timestepper_detail::assemble_helper<ExecSpace, DerivType>::assemble_k_total(
                    k_total,
                    KOKKOS_LAMBDA(std::array<DerivType, 2> k) { return k[0] + k[1]; },
                    k1,
                    k_new);

            // Save the old characteristic feet
            timestepper_detail::copy_helper<ValType>::copy(y_old, y);

            // Re-initialise the characteristic feet
            timestepper_detail::copy_helper<ValType>::copy(y, y_init);

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, k_total, 0.5 * dt);

            // Check convergence
            not_converged = (norm_inf(y_old - y) / norm_inf(y_old)) >= m_epsilon;

        } while (not_converged and (counter < m_max_counter));
    }
};

template <
        timestepper_detail::FieldLike FieldMem,
        timestepper_detail::FieldLike DerivFieldMem,
        class ExecSpace>
class CrankNicolson<FieldMem, DerivFieldMem, ExecSpace>
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
        using element_type = typename timestepper_detail::ElementType<DerivField>::type;

        FieldMem y_init_alloc("y_init (CrankNicholson::update)", m_idx_range);
        FieldMem y_old_alloc("y_old (CrankNicholson::update)", m_idx_range);
        DerivFieldMem k1_alloc("k1 (CrankNicholson::update)", m_idx_range);
        DerivFieldMem k_new_alloc("k_new (CrankNicholson::update)", m_idx_range);
        DerivFieldMem k_total_alloc("k_total (CrankNicholson::update)", m_idx_range);

        ValField y_init = get_field(y_init_alloc);
        ValField y_old = get_field(y_old_alloc);
        DerivField k1 = get_field(k1_alloc);
        DerivField k_new = get_field(k_new_alloc);
        DerivField k_total = get_field(k_total_alloc);

        // Save initial conditions
        timestepper_detail::copy_helper<FieldMem>::copy(y_init, get_const_field(y));

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
            timestepper_detail::assemble_helper<ExecSpace, DerivFieldMem>::assemble_k_total(
                    exec_space,
                    k_total,
                    KOKKOS_LAMBDA(std::array<element_type, 2> k) { return k[0] + k[1]; },
                    k1,
                    k_new);

            // Save the old characteristic feet
            timestepper_detail::copy_helper<FieldMem>::copy(y_old, get_const_field(y));

            // Re-initialise the characteristic feet
            timestepper_detail::copy_helper<FieldMem>::copy(y, get_const_field(y_init));

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, get_const_field(k_total), 0.5 * dt);


            // Check convergence
            not_converged = (error_norm_inf(exec_space, get_const_field(y_old), get_const_field(y))
                             / norm_inf(exec_space, get_const_field(y_old)))
                            >= m_epsilon;


        } while (not_converged and (counter < m_max_counter));
    }
};

class CrankNicolsonBuilder
{
private:
    int const m_max_counter;
    double const m_epsilon;

public:
    explicit CrankNicolsonBuilder(int const counter = 20, double const epsilon = 1e-12)
        : m_max_counter(counter)
        , m_epsilon(epsilon)
    {
    }

    template <
            class FieldMem,
            class DerivFieldMem = FieldMem,
            class ExecSpace = Kokkos::DefaultExecutionSpace>
    using time_stepper_t = CrankNicolson<FieldMem, DerivFieldMem, ExecSpace>;

    template <class TimeStepper>
    auto preallocate(typename TimeStepper::IdxRange const idx_range) const
    {
        static_assert(std::is_same_v<
                      TimeStepper,
                      time_stepper_t<
                              typename TimeStepper::ValFieldMem,
                              typename TimeStepper::DerivFieldMem,
                              typename TimeStepper::exec_space>>);
        return TimeStepper(idx_range, m_max_counter, m_epsilon);
    }
};

namespace detail {

template <>
inline constexpr bool enable_is_timestepper_builder<CrankNicolsonBuilder> = true;

}
```


