// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <type_traits>

#include <ddc_helper.hpp>
#include <vector_field_common.hpp>

#include "ddc_aliases.hpp"
#include "itimestepper.hpp"

/**
 * @brief A class which provides an implementation of a second-order Runge-Kutta method.
 *
 * A class which provides an implementation of a second-order Runge-Kutta method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
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
template <class FieldMemType, class DerivFieldMemType = FieldMemType>
class RK2 : public ITimeStepper
{
private:
    static_assert(ddc::is_chunk_v<FieldMemType> or is_field_v<FieldMemType>);
    static_assert(ddc::is_chunk_v<DerivFieldMemType> or is_field_v<DerivFieldMemType>);

    static_assert(std::is_same_v<
                  typename FieldMemType::discrete_domain_type,
                  typename DerivFieldMemType::discrete_domain_type>);

    using IdxRange = typename FieldMemType::discrete_domain_type;

    using Idx = typename IdxRange::discrete_element_type;

    using ValField = typename FieldMemType::span_type;
    using ValConstField = typename FieldMemType::view_type;

    using DerivField = typename DerivFieldMemType::span_type;
    using DerivConstField = typename DerivFieldMemType::view_type;


    IdxRange const m_dom;

public:
    /**
     * @brief Create a RK2 object.
     * @param[in] dom The index range on which the points which evolve over time are defined.
     */
    explicit RK2(IdxRange dom) : m_dom(dom) {}

    /**
     * @brief Carry out one step of the Runge-Kutta scheme.
     *
     * This function is a wrapper around the update function below. The values of the function are
     * updated using the trivial method $f += df * dt$. This is the standard method however some
     * cases may need a more complex update function which is why the more explicit method is
     * also provided.
     *
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the index range.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy
     *     The function describing how the derivative of the evolve function is calculated.
     */
    void update(ValField y, double dt, std::function<void(DerivField, ValConstField)> dy) const
    {
        using ExecSpace = typename FieldMemType::memory_space::execution_space;
        update(ExecSpace(), y, dt, dy);
    }

    /**
     * @brief Carry out one step of the Runge-Kutta scheme.
     *
     * This function is a wrapper around the update function below. The values of the function are
     * updated using the trivial method $f += df * dt$. This is the standard method however some
     * cases may need a more complex update function which is why the more explicit method is
     * also provided.
     *
     * @param[in] exec_space
     *     The space on which the function is executed (CPU/GPU).
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the index range.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy
     *     The function describing how the derivative of the evolve function is calculated.
     */
    template <class ExecSpace>
    void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy) const
    {
        static_assert(ddc::is_chunk_v<FieldMemType>);
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename FieldMemType::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename DerivFieldMemType::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        update(exec_space, y, dt, dy, [&](ValField y, DerivConstField dy, double dt) {
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(y),
                    KOKKOS_LAMBDA(Idx const idx) { y(idx) = y(idx) + dy(idx) * dt; });
        });
    }

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
     * @param[in] dy
     *     The function describing how the derivative of the evolve function is calculated.
     * @param[in] y_update
     *     The function describing how the value(s) are updated using the derivative.
     */
    template <class ExecSpace>
    void update(
            ExecSpace const& exec_space,
            ValField y,
            double dt,
            std::function<void(DerivField, ValConstField)> dy,
            std::function<void(ValField, DerivConstField, double)> y_update) const
    {
        DerivFieldMemType m_k1(m_dom);
        DerivFieldMemType m_k2(m_dom);
        FieldMemType m_y_prime(m_dom);


        // Save initial conditions
        if constexpr (is_field_v<FieldMemType>) {
            ddcHelper::deepcopy(m_y_prime, y);
        } else {
            ddc::parallel_deepcopy(m_y_prime, y);
        }

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y)
        dy(m_k1, y);

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(m_y_prime, m_k1, 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy(m_k2, m_y_prime);

        // ----------- Update y --------------
        // Calculate y_{n+1} := y_n + h*k_2
        y_update(y, m_k2, dt);
    }
};
