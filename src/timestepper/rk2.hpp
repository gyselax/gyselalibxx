// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "vector_field_common.hpp"

/**
 * @brief A class which provides an implementation of a second-order Runge-Kutta method.
 *
 * A class which provides an implementation of a second-order Runge-Kutta method in
 * order to evolve values over time. This specialisation handles elementwise
 * operations and can be called from GPU.
 *
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the Runge-Kutta 2 method is given by :
 * @f$ y^{n+1} =  y^{n} + dt k_2 @f$,
 *
 * with
 *
 * - @f$ k_1 = f(t^{n}, y^{n}) @f$,
 * - @f$ k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 ) @f$,
 */
template <class ValType, class DerivType = ValType, class ExecSpace = Kokkos::DefaultExecutionSpace>
class RK2
{
    static_assert(!timestepper_detail::FieldLike<ValType>);
    static_assert(!timestepper_detail::FieldLike<DerivType>);

public:
    /**
     * @brief Create a RK2 object to operate on scalars.
     */
    explicit KOKKOS_DEFAULTED_FUNCTION RK2() = default;

    /**
     * @brief Carry out one step of the Runge-Kutta scheme on a scalar.
     *
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
    template <
            class DYFunctor,
            class YFunctor
            = decltype(timestepper_detail::serial_y_updater<ValType&, DerivType const&>::y_update)>
    KOKKOS_FUNCTION void update(
            ValType y,
            double dt,
            DYFunctor dy_calculator,
            YFunctor y_update
            = timestepper_detail::serial_y_updater<ValType&, DerivType const&>::y_update) const
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

/**
 * @brief A class which provides an implementation of a second-order Runge-Kutta method.
 *
 * A class which provides an implementation of a second-order Runge-Kutta method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * This specialisation handles Field-like objects (Field, VectorField, MultipatchField).
 * The values which evolve are defined on an index range.
 *
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the Runge-Kutta 2 method is given by :
 * @f$ y^{n+1} =  y^{n} + dt k_2 @f$,
 *
 * with
 *
 * - @f$ k_1 = f(t^{n}, y^{n}) @f$,
 * - @f$ k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 ) @f$,
 */
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
    /**
     * @brief Create a RK2 object.
     * @param[in] idx_range The index range on which the points which evolve over time are defined.
     */
    explicit RK2(IdxRange idx_range) : m_idx_range(idx_range) {}

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
