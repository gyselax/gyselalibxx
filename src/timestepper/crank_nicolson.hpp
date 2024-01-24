#pragma once
#include <ddc_helper.hpp>
#include <vector_field_common.hpp>

#include "itimestepper.hpp"
#include "utils_tools.hpp"


/**
 * @brief A class which provides an implementation of a Crank-Nicolson method.
 *
 * A class which provides an implementation of a Crank-Nicolson method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * The values which evolve are defined on a domain.
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
template <class ValChunk, class DerivChunk = ValChunk>
class CrankNicolson : public ITimeStepper
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
    int const m_max_counter;
    double const m_epsilon;

public:
    /**
     * @brief Create a CrankNicolson object.
     * @param[in] dom
     *      The domain on which the points which evolve over time are defined.
     * @param[in] counter
     *      The maximal number of loops for the implicit method.
     * @param[in] epsilon
     *      The @f$ \varepsilon @f$ upperbound of the difference of two steps
     *      in the implicit method: @f$ |y^{k+1} -  y^{k}| < \varepsilon @f$.
     */
    CrankNicolson(Domain dom, int const counter = int(20), double const epsilon = 1e-12)
        : m_dom(dom)
        , m_max_counter(counter)
        , m_epsilon(epsilon)
    {
    }

    /**
     * @brief Carry out one step of the Crank-Nicolson scheme.
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
     * @brief Carry out one step of the Crank-Nicolson scheme.
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
        static_assert(ddc::is_chunk_v<ValChunk>);
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename ValChunk::memory_space>::accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        static_assert(
                Kokkos::SpaceAccessibility<ExecSpace, typename DerivChunk::memory_space>::
                        accessible,
                "MemorySpace has to be accessible for ExecutionSpace.");
        update(exec_space, y, dt, dy, [&](ValSpan y, DerivView dy, double dt) {
            ddc::for_each(
                    ddc::policies::policy(exec_space),
                    y.domain(),
                    KOKKOS_LAMBDA(Index const idx) { y(idx) = y(idx) + dy(idx) * dt; });
        });
    }

    /**
     * @brief Carry out one step of the Crank-Nicolson scheme.
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
        ValChunk m_y_init_alloc(m_dom);
        ValChunk m_y_old_alloc(m_dom);
        DerivChunk m_k1_alloc(m_dom);
        DerivChunk m_k_new_alloc(m_dom);
        DerivChunk m_k_total_alloc(m_dom);
        ValSpan m_y_init = m_y_init_alloc.span_view();
        ValSpan m_y_old = m_y_old_alloc.span_view();
        DerivSpan m_k1 = m_k1_alloc.span_view();
        DerivSpan m_k_new = m_k_new_alloc.span_view();
        DerivSpan m_k_total = m_k_total_alloc.span_view();


        copy(m_y_init, y);

        // --------- Calculate k1 ------------
        // Calculate k1 = f(y_n)
        dy(m_k1, y);

        // -------- Calculate k_new ----------
        bool not_converged = true;
        int counter = 0;
        do {
            counter++;

            // Calculate k_new = f(y_new)
            dy(m_k_new, y);

            // Calculation of step
            if constexpr (is_field_v<DerivChunk>) {
                ddc::for_each(
                        ddc::policies::policy(exec_space),
                        m_k_total.domain(),
                        KOKKOS_CLASS_LAMBDA(Index const i) {
                            // k_total = k1 + k_new
                            fill_k_total(i, m_k_total, m_k1(i) + m_k_new(i));
                        });
            } else {
                ddc::for_each(
                        ddc::policies::policy(exec_space),
                        m_k_total.domain(),
                        KOKKOS_CLASS_LAMBDA(Index const i) {
                            // k_total = k1 + k_new
                            m_k_total(i) = m_k1(i) + m_k_new(i);
                        });
            }

            // Save the old characteristic feet
            copy(m_y_old, y);

            // Re-initiliase the characteristic feet
            copy(y, m_y_init);

            // Calculate y_new := y_n + h/2*(k_1 + k_new)
            y_update(y, m_k_total, 0.5 * dt);


            // Check convergence
            not_converged = not have_converged(exec_space, m_y_old, y);


        } while (not_converged and (counter < m_max_counter));
    };


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
    template <class ExecSpace>
    bool have_converged(ExecSpace const& exec_space, ValView y_old, ValView y_new) const
    {
        auto const dom = y_old.domain();

        double norm_old = ddc::transform_reduce(
                ddc::policies::policy(exec_space),
                dom,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(Index const idx) { return norm_inf(y_old(idx)); });

        double max_diff = ddc::transform_reduce(
                ddc::policies::policy(exec_space),
                dom,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(Index const idx) { return norm_inf(y_old(idx) - y_new(idx)); });

        return (max_diff / norm_old) < m_epsilon;
    }

private:
    void copy(ValSpan copy_to, ValView copy_from) const
    {
        if constexpr (is_field_v<ValSpan>) {
            ddcHelper::deepcopy(copy_to, copy_from);
        } else {
            ddc::deepcopy(copy_to, copy_from);
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
