#pragma once
#include <array>

#include <ddc_helper.hpp>
#include <vector_field_common.hpp>

#include "itimestepper.hpp"

/**
 * @brief A class which provides an implementation of a second-order Runge-Kutta method.
 *
 * A class which provides an implementation of a second-order Runge-Kutta method in
 * order to evolve values over time. The values may be either scalars or vectors. In the
 * case of vectors the appropriate dimensions must be passed as template parameters.
 * The values which evolve are defined on a domain.
 */
template <class ValChunk, class DerivChunk = ValChunk>
class RK2 : public ITimeStepper
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

    DerivChunk m_k1;
    DerivChunk m_k2;
    ValChunk m_y_prime;

public:
    /**
     * @brief Create a RK2 object.
     * @param[in] dom The domain on which the points which evolve over time are defined.
     */
    RK2(Domain dom)
    {
        m_k1 = DerivChunk(dom);
        m_k2 = DerivChunk(dom);
        m_y_prime = ValChunk(dom);
    }

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
    void update(ValSpan y, double dt, std::function<void(DerivSpan, ValView)> dy)
    {
        static_assert(ddc::is_chunk_v<ValChunk>);
        update(y, dt, dy, [&](ValSpan y, DerivView dy, double dt) {
            ddc::for_each(y.domain(), [&](Index const idx) { y(idx) = y(idx) + dy(idx) * dt; });
        });
    }

    /**
     * @brief Carry out one step of the Runge-Kutta scheme.
     *
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
    void update(
            ValSpan y,
            double dt,
            std::function<void(DerivSpan, ValView)> dy,
            std::function<void(ValSpan, DerivView, double)> y_update)
    {
        // Save initial conditions
        if constexpr (is_field_v<ValChunk>) {
            ddcHelper::deepcopy(m_y_prime, y);
        } else {
            ddc::deepcopy(m_y_prime, y);
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
    };
};
