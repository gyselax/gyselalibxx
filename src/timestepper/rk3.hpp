#pragma once
#include <ddc_helper.hpp>
#include <vector_field_common.hpp>

#include "itimestepper.hpp"


/**
 * @brief A class which provides an implementation of a third-order Runge-Kutta method.
 *
 * A class which provides an implementation of a third-order Runge-Kutta method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * The values which evolve are defined on a domain.
 *
 * For the following ODE :
 * @f$\partial_t y(t) = f(t, y(t)) @f$,
 *
 * the Runge-Kutta 3 method is given by :
 * @f$ y^{n+1} =  y^{n} + \frac{dt}{6} \left(k_1 + 4 k_2 + k_3 \right) @f$,
 *
 * with
 *
 * - @f$ k_1 = f(t^{n}, y^{n}) @f$,
 * - @f$ k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 ) @f$,
 * - @f$ k_3 = f(t^{n+1/2}, y^{n} + dt ( 2 k_2 - k_1) ) @f$.
 *
 */
template <class ValChunk, class DerivChunk = ValChunk>
class RK3 : public ITimeStepper
{
private:
    static_assert(ddc::is_chunk_v<ValChunk> or is_field_v<ValChunk>);
    static_assert(ddc::is_chunk_v<DerivChunk> or is_field_v<DerivChunk>);

    static_assert(
            std::is_same_v<typename ValChunk::mdomain_type, typename DerivChunk::mdomain_type>);

    using Domain = typename ValChunk::mdomain_type;

    using Index = typename Domain::discrete_element_type;

    using ValSpan = typename ValChunk::span_type;
    using ValView = typename ValChunk::view_type;

    using DerivSpan = typename DerivChunk::span_type;
    using DerivView = typename DerivChunk::view_type;

    Domain const m_dom;

public:
    /**
     * @brief Create a RK3 object.
     * @param[in] dom The domain on which the points which evolve over time are defined.
     */
    RK3(Domain dom) : m_dom(dom) {}

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
     *     of the domain.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy
     *     The function describing how the derivative of the evolve function is calculated.
     */
    void update(ValSpan y, double dt, std::function<void(DerivSpan, ValView)> dy) const
    {
        using ExecSpace = typename ValChunk::memory_space::execution_space;
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
     *     of the domain.
     * @param[in] dt
     *     The time step over which the values should be evolved.
     * @param[in] dy
     *     The function describing how the derivative of the evolve function is calculated.
     */
    template <class ExecSpace>
    void update(
            ExecSpace const& exec_space,
            ValSpan y,
            double dt,
            std::function<void(DerivSpan, ValView)> dy) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename ValChunk::memory_space>::accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename DerivChunk::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(ddc::is_chunk_v<ValChunk>);
        update(exec_space, y, dt, dy, [&](ValSpan y, DerivView dy, double dt) {
            ddc::parallel_for_each(
                    exec_space,
                    y.domain(),
                    KOKKOS_LAMBDA(Index const idx) { y(idx) = y(idx) + dy(idx) * dt; });
        });
    }

    /**
     * @brief Carry out one step of the Runge-Kutta scheme.
     *
     * @param[in] exec_space
     *     The space on which the function is executed (CPU/GPU).
     * @param[inout] y
     *     The value(s) which should be evolved over time defined on each of the dimensions at each point
     *     of the domain.
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
            ValSpan y,
            double dt,
            std::function<void(DerivSpan, ValView)> dy,
            std::function<void(ValSpan, DerivView, double)> y_update) const
    {
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename ValChunk::memory_space>::accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename DerivChunk::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");

        ValChunk m_y_prime_alloc(m_dom);
        DerivChunk m_k1_alloc(m_dom);
        DerivChunk m_k2_alloc(m_dom);
        DerivChunk m_k3_alloc(m_dom);
        DerivChunk m_k_total_alloc(m_dom);
        ValSpan m_y_prime = m_y_prime_alloc.span_view();
        DerivSpan m_k1 = m_k1_alloc.span_view();
        DerivSpan m_k2 = m_k2_alloc.span_view();
        DerivSpan m_k3 = m_k3_alloc.span_view();
        DerivSpan m_k_total = m_k_total_alloc.span_view();

        // Save initial conditions
        copy(m_y_prime, y);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y)
        dy(m_k1, y);

        // --------- Calculate k2 ------------
        // Calculate y_new := y_n + h/2*k_1
        y_update(m_y_prime, m_k1, 0.5 * dt);

        // Calculate k2 = f(y_new)
        dy(m_k2, m_y_prime);

        // --------- Calculate k3 ------------
        // Calculation of step
        if constexpr (is_field_v<DerivChunk>) {
            ddc::parallel_for_each(
                    exec_space,
                    m_k_total.domain(),
                    KOKKOS_CLASS_LAMBDA(Index const i) {
                        // k_total = 2 * k2 - k1
                        fill_k_total(i, m_k_total, 2 * m_k2(i) - m_k1(i));
                    });
        } else {
            ddc::parallel_for_each(
                    exec_space,
                    m_k_total.domain(),
                    KOKKOS_LAMBDA(Index const i) {
                        // k_total = 2 * k2 - k1
                        m_k_total(i) = 2 * m_k2(i) - m_k1(i);
                    });
        }

        // Collect initial conditions
        copy(m_y_prime, y);

        // Calculate y_new := y_n + h*(2*k_2-k_1)
        y_update(m_y_prime, m_k_total, dt);

        // Calculate k3 = f(y_new)
        dy(m_k3, m_y_prime);

        // --------- Update y ------------
        // Calculation of step
        if constexpr (is_field_v<DerivChunk>) {
            ddc::parallel_for_each(
                    exec_space,
                    m_k_total.domain(),
                    KOKKOS_CLASS_LAMBDA(Index const i) {
                        // k_total = k1 + 4 * k2 + k3
                        fill_k_total(i, m_k_total, m_k1(i) + 4 * m_k2(i) + m_k3(i));
                    });
        } else {
            ddc::parallel_for_each(
                    exec_space,
                    m_k_total.domain(),
                    KOKKOS_LAMBDA(Index const i) {
                        // k_total = k1 + 4 * k2 + k3
                        m_k_total(i) = m_k1(i) + 4 * m_k2(i) + m_k3(i);
                    });
        }

        // Calculate y_{n+1} := y_n + (k1 + 4 * k2 + k3) * h/6
        y_update(y, m_k_total, dt / 6.);
    }

private:
    void copy(ValSpan copy_to, ValView copy_from) const
    {
        if constexpr (is_field_v<ValSpan>) {
            ddcHelper::deepcopy(copy_to, copy_from);
        } else {
            ddc::parallel_deepcopy(copy_to, copy_from);
        }
    }

    template <class... DDims>
    KOKKOS_FUNCTION void fill_k_total(
            Index i,
            DerivSpan m_k_total,
            ddc::Coordinate<DDims...> new_val) const
    {
        ((ddcHelper::get<DDims>(m_k_total)(i) = ddc::get<DDims>(new_val)), ...);
    }
};
