// SPDX-License-Identifier: MIT
#pragma once
#include <array>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "vector_field_common.hpp"

/**
 * @brief A class which provides an implementation of an explicit Euler method.
 *
 * A class which provides an implementation of an explicit Euler method in
 * order to evolve values over time. This specialisation handles elementwise
 * operations and can be called from GPU.
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the explicit Euler method is given by :
 * @f$ y^{n+1} =  y^{n} + dt f(t^{n}, y^{n})@f$.
 *
 * The method is order 1.
 *
 */
template <class ValType, class DerivType = ValType, class ExecSpace = Kokkos::DefaultExecutionSpace>
class Euler
{
    static_assert(!timestepper_detail::FieldLike<ValType>);
    static_assert(!timestepper_detail::FieldLike<DerivType>);

public:
    /// The type of the memory allocation for the values of the function being evolved.
    using ValFieldMem = ValType;

    /// The type of the memory allocation for the derivatives of the function being evolved.
    using DerivFieldMem = DerivType;

    /// The space (CPU/GPU) where the calculations are carried out.
    using exec_space = ExecSpace;

public:
    /**
     * @brief Create a Euler object to operate on scalars.
     */
    explicit KOKKOS_DEFAULTED_FUNCTION Euler() = default;

    /**
     * @brief Carry out one step of the explicit Euler scheme on a scalar.
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
            = decltype(timestepper_detail::serial_y_update<ValType&, DerivType const&>)>
    KOKKOS_FUNCTION void update(
            ValType y,
            double dt,
            DYFunctor dy_calculator,
            YFunctor y_update
            = timestepper_detail::serial_y_update<ValType&, DerivType const&>) const
    {
        static_assert(std::is_invocable_v<DYFunctor, DerivType&, ValType>);
        DerivType k1;

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(k1, y);

        // ----------- Update y --------------
        // Calculate y_new := y_n + h*k_1
        y_update(y, k1, dt);
    }
};

/**
 * @brief A class which provides an implementation of an explicit Euler method.
 *
 * A class which provides an implementation of an explicit Euler method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * This specialisation handles Field-like objects (Field, VectorField, MultipatchField).
 * The values which evolve are defined on an index range.
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the explicit Euler method is given by :
 * @f$ y^{n+1} =  y^{n} + dt f(t^{n}, y^{n})@f$.
 *
 * The method is order 1.
 *
 */
template <
        timestepper_detail::FieldLike FieldMem,
        timestepper_detail::FieldLike DerivFieldMem,
        class ExecSpace>
class Euler<FieldMem, DerivFieldMem, ExecSpace>
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
     * @brief Create a Euler object.
     * @param[in] idx_range The index range on which the points which evolve over time are defined.
     */
    explicit Euler(IdxRange idx_range) : m_idx_range(idx_range) {}

    /**
     * @brief Carry out one step of the explicit Euler scheme.
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
        DerivFieldMem k1_alloc("k1 (Euler::update)", m_idx_range);
        DerivField k1 = get_field(k1_alloc);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(k1, get_const_field(y));

        // ----------- Update y --------------
        // Calculate y_new := y_n + h*k_1
        y_update(y, get_const_field(k1), dt);
    }
};

using EulerBuilder = ExplicitTimeStepperBuilder<Euler>;
