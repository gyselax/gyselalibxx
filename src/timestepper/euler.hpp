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
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
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
        class FieldMem,
        class DerivFieldMem = FieldMem,
        class ExecSpace = Kokkos::DefaultExecutionSpace>
class Euler : public ITimeStepper<FieldMem, DerivFieldMem, ExecSpace>
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
        DerivFieldMem k1_alloc(m_idx_range);
        DerivField k1(k1_alloc);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy_calculator(k1, get_const_field(y));

        // ----------- Update y --------------
        // Calculate y_new := y_n + h*k_1
        y_update(y, get_const_field(k1), dt);
    }
};
