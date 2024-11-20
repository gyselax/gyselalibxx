// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief The diocotron exact solution of the density @f$ \rho @f$.
 *
 * The equations of the diocotron simulation are the Vlasov-Poisson equations
 *
 * - @f$ \partial_t \rho - E_y \partial_x \rho + E_x \partial_y \rho = 0 @f$,
 * - @f$ L \phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho @f$,
 * - @f$ E = -\nabla \phi @f$.
 *
 * The DiocotronDensitySolution provides an initial perturbed solution @f$ \rho_0 @f$ for
 * the density @f$ \rho @f$ and its associated equilibrium solution @f$ \rho_{eq} @f$.
 */
class DiocotronDensitySolution
{
private:
    /**
     * @brief The value of @f$ r @f$ at the inner wall.
     */
    double const m_W1;
    /**
     * @brief The value of @f$ r @f$ at the inner boundary
     * of the initial condition
     */
    double const m_R1;
    /**
     * @brief The value of @f$ r @f$ at the outer boundary
     * of the initial condition
     */
    double const m_R2;
    /**
     * @brief The value of @f$ r @f$ at the outer wall.
     */
    double const m_W2;

    /**
     * @brief The charge carried by the inner conductor wall.
     */
    double const m_Q;
    /**
     * @brief The mode of the perturbation.
     */
    int const m_l;
    /**
     * @brief The amplitude of the perturbation.
     */
    double const m_eps;

    /**
     * @brief The real part of @f$ \omega @f$,
     * solution of the dispersion relation.
     */
    double m_omega_Re;
    /**
     * @brief The imaginary part of @f$ \omega @f$,
     * solution of the dispersion relation.
     */
    double m_omega_Im;

    /**
     * @brief The average between @f$ R_1 @f$ and @f$ R_2 @f$.
     */
    const double m_r_bar;
    /**
     * @brief The half difference between @f$ R_1 @f$ and @f$ R_2 @f$.
     */
    const double m_d;
    /**
     * @brief The power inside the exponential function.
     */
    const double m_p;


public:
    /**
     * @brief Instantiate a DiocotronDensitySolution.
     *
     *
     * @param[in] W1
     *      The value of the r-coordinate where the index range starts.
     * @param[in] R1
     *      The value of the r-coordinate where the index range of non null
     *      initial condition starts.
     * @param[in] R2
     *      The value of the r-coordinate where the index range of non null
     *      initial condition ends.
     * @param[in] W2
     *      The value of the r-coordinate where the index range ends.
     * @param[in] Q
     *      The charge carried by the inner conductor at @f$ r = W_1 @f$.
     * @param[in] l
     *      The mode of the perturbation.
     * @param[in] eps
     *      The amplitude of the perturbation.
     */
    DiocotronDensitySolution(
            double const W1,
            double const R1,
            double const R2,
            double const W2,
            double const Q,
            int const l,
            double const eps)
        : m_W1(W1)
        , m_R1(R1)
        , m_R2(R2)
        , m_W2(W2)
        , m_Q(Q)
        , m_l(l)
        , m_eps(eps)
        , m_omega_Re(0.)
        , m_omega_Im(0.)
        , m_r_bar((R1 + R2) / 2.)
        , m_d((R2 - R1) / 2.)
        , m_p(50.)
    {
        assert(0 <= m_W1);
        assert(m_W1 <= m_R1);
        assert(m_R1 < m_R2);
        assert(m_R2 <= m_W2);
        assert(m_l != 0);

        // Computation of omega: -----------------------------------------------------------------
        dispersion_relation(m_omega_Re, m_omega_Im);
    };

    ~DiocotronDensitySolution() {};

    /**
     * @brief Get the frequency of the perturbation.
     *
     * The frequency of the perturbation is given by the
     * real part of @f$ \omega @f$, where @f$ \omega @f$
     * is the solution of the dispersion relation.
     *
     * @return The frequency of the perturbation.
     */
    double get_frequency() const;

    /**
     * @brief Get the slope of the perturbation.
     *
     * The slope of the perturbation is given by the
     * imaginary part of @f$ \omega @f$, where @f$ \omega @f$
     * is the solution of the dispersion relation.
     *
     * @return The slope of the perturbation.
     */
    double get_slope() const;


    /**
     * @brief Get the initial condition of the density.
     *
     * The initial condition is given by
     *
     * - if @f$ R_1 \leq r \leq R_2 @f$, @f$ \rho_0(r, \theta) = (1 + \varepsilon \cos(l\theta)) \exp\left( - \left(\frac{r - \bar{r}}{d}\right)^p \right) @f$,
     * - otherwise, @f$ \rho_0(r, \theta) = 0 @f$,
     *
     * with @f$ p = 50 @f$, @f$ \bar{r} = \frac{R_1 + R_2}{2}@f$ and  @f$ d = \frac{R_2 - R_1}{2}@f$.
     *
     * @param[in] coord
     *      The coordinate where we evaluate the initial condition.
     *
     * @return The value of the initial condition at the given coordinate.
     */
    double initialisation(Coord<R, Theta> const& coord) const;

    /**
     * @brief Get the equilibrium of the density.
     *
     * The equilibrium is given by
     *
     * - if @f$ R_1 \leq r \leq R_2 @f$, @f$ \rho_0(r, \theta) = \exp\left( - \left(\frac{r - \bar{r}}{d}\right)^p \right) @f$,
     * - otherwise, @f$ \rho_0(r, \theta) = 0 @f$,
     *
     * with @f$ p = 50 @f$, @f$ \bar{r} = \frac{R_1 + R_2}{2}@f$ and @f$ d = \frac{R_2 - R_1}{2}@f$.
     *
     * @param[in] coord
     *      The coordinate where we evaluate the equilibrium.
     *
     * @return The value of the equilibrium at the given coordinate.
     */
    double equilibrium(Coord<R, Theta> const& coord) const;

private:
    /**
     * @brief Solves the dispersion relation.
     *
     * The dispersion relation is given by
     * @f$ \left(\frac{\omega}{\omega_D} \right)^2 - b_l \frac{\omega}{\omega_D} + c_l = 0 @f$
     *
     * with
     * - @f$ b_l = \frac{1}{1 - \left(\frac{W_1}{W_2}\right)^{2l}}
     *      \left(
     *          l \left[ 1 - \left(\frac{W_1}{W_2}\right)^{2l} \right]
     *          \left(
     *              \left[ 1 - \left(\frac{R_1}{R_2}\right)^{2} \right]
     *              + \frac{\omega_q}{\omega_D}\left[1 +  \left(\frac{R_1}{R_2}\right)^{2} \right]
     *          \right)
     *          +
     *           \left[ 1 - \left(\frac{R_1}{R_2}\right)^{2l} \right]
     *           \left[\left(\frac{R_2}{W_2}\right)^{2l} - \left(\frac{W_1}{R_1}\right)^{2l} \right]
     *      \right)@f$,
     *
     * - @f$ c_l = \frac{1}{1 - \left(\frac{W_1}{W_2}\right)^{2l}}
     *      \left(
     *          l^2 \frac{\omega_q}{\omega_D} \left[1 - \left(\frac{R_1}{R_2}\right)^{2} + \frac{\omega_q}{\omega_D}\left(\frac{R_1}{R_2}\right)^{2} \right]
     *              \left[ 1 - \left(\frac{W_1}{W_2}\right)^{2l} \right]
     *          - l \frac{\omega_q}{\omega_D} \left[ 1 - \left(\frac{W_1}{R_2}\right)^{2l} \right]
     *              \left[ 1 - \left(\frac{R_2}{W_2}\right)^{2l} \right]
     *          + l  \left[1 - \left(\frac{R_1}{R_2}\right)^{2} + \frac{\omega_q}{\omega_D}\left(\frac{R_1}{R_2}\right)^{2} \right]
     *              \left[ 1 - \left(\frac{W_1}{R_1}\right)^{2l} \right] \left[ 1 - \left(\frac{R_1}{W_2}\right)^{2l} \right]
     *          -  \left[ 1 - \left(\frac{R_2}{W_2}\right)^{2l} \right]  \left[ 1 - \left(\frac{W_1}{R_1}\right)^{2l} \right]
     *               \left[ 1 - \left(\frac{R_1}{R_2}\right)^{2l} \right]
     *      \right)@f$,
     *
     * and
     * - @f$ \omega_D = 0.5 @f$,
     * - @f$ \omega_q = - \frac{2 Q c}{B_0 R_1^2}@f$ (here @f$ c = 1 @f$, @f$ B_0 = 1 @f$).
     *
     * (For more details, see R.C. Davidson's book, "Theory of Nonneutral Plasmas",
     * Chap.6, "The diocotron instability", Benjamin, 1974)
     *
     * @param[out] omega_Re
     *      The real part of @f$ \omega @f$.
     * @param[out] omega_Im
     *      The imaginary part of @f$ \omega @f$.
     */
    void dispersion_relation(double& omega_Re, double& omega_Im) const;
};
