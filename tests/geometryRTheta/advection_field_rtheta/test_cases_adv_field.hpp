// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "../../advection/r_theta_test_cases.hpp"

#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "math_tools.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"

/*
 *  This file defines
 *   - the TEST ELECTROSTATIC POTENTIALS,
 *   - and the TEST SIMULATIONS.
 */



// TEST ELECTROSTATIC POTENTIALS -----------------------------------------------------------------
/**
 * @brief Electrostatic potential for a decentred rotation test of the 2D polar advection operator.
 *
 * The test advection field is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889):
 *
 * @f$ A(x,y) = \omega [y_c - y, x - x_c]^T = [-\partial_y \phi, \partial_x \phi]^T @f$,
 * 
 * An adapted electrostatic potential is 
 * @f$ \phi(x,y) = \omega (1\2 * x^2 - x x_c + 1/2 * y^2 - y y_c)@f$.
 *
 * The characteristic feet are then given by
 * - @f$ X(t + dt) = x_c + (X(t) - x_c) \cos(\omega dt) - (Y(t) - y_c) \sin(\omega dt) @f$,
 * - @f$ Y(t + dt) = y_c + (X(t) - x_c) \sin(\omega dt) + (Y(t) - y_c) \cos(\omega dt) @f$.
 *
 */
class ElectrostaticPotentialSimulation_decentred_rotation
{
private:
    double const m_omega;
    double const m_xc;
    double const m_yc;

public:
    /// @brief Instantiate an ElectrostaticPotentialSimulation_decentred_rotation advection field.
    KOKKOS_FUNCTION ElectrostaticPotentialSimulation_decentred_rotation()
        : m_omega(2 * M_PI)
        , m_xc(0.25)
        , m_yc(0.)
    {
    }

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION ElectrostaticPotentialSimulation_decentred_rotation(
            ElectrostaticPotentialSimulation_decentred_rotation const&)
            = default;

    /**
     * @brief Get the electrostatic potential associated to the advection field.
     *
     * @param[in] coord
     *      The coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The electrostatic potential at the given coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordXY const coord, double const t) const
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        return m_omega * (0.5 * x * x - m_xc * x + 0.5 * y * y - m_yc * y);
    }

    /**
     * @brief Get the characteristic feet in the physical domain.
     *
     * @param[in] coord
     *      The original coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The characteristic feet in the physical domain.
     */
    KOKKOS_FUNCTION CoordXY exact_feet(CoordXY coord, double const t) const
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        double const foot_x = m_xc + (x - m_xc) * Kokkos::cos(m_omega * -t)
                              - (y - m_yc) * Kokkos::sin(m_omega * -t);
        double const foot_y = m_yc + (x - m_xc) * Kokkos::sin(m_omega * -t)
                              + (y - m_yc) * Kokkos::cos(m_omega * -t);
        return CoordXY(foot_x, foot_y);
    }
};



/**
 * @brief Electrostatic potential for a translation test of the 2D polar advection operator.
 *
 * The test advection field for a translation in the physical domain is given by :
 *
 * @f$ A(x,y) = [v_x, v_y] = [-\partial_y \phi, \partial_x \phi]^T @f$,
 *
 * with @f$ v_x @f$ and @f$ v_y @f$ constants.
 *
 * An adapted electrostatic potential is 
 * @f$ \phi(x,y) = v_y x - v_x y@f$.
 * 
 * The characteristic feet are then given by
 * - @f$ X(t + dt) = X(t) + dt v_x @f$,
 * - @f$ Y(t + dt) = Y(t) + dt v_y @f$.
 *
 */
class ElectrostaticPotentialSimulation_translation
{
private:
    CoordXY const m_velocity;

public:
    /**
     * @brief Instantiate an ElectrostaticPotentialSimulation_translation advection field.
     *
     * @param[in] vx
     *      The constant first component of the advection field in the physical domain.
     * @param[in] vy
     *      The constant second component of the advection field in the physical domain.
     */
    ElectrostaticPotentialSimulation_translation(CoordVx vx, CoordVy vy)
        : m_velocity(ddc::get<Vx>(vx), ddc::get<Vy>(vy))
    {
    }

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION ElectrostaticPotentialSimulation_translation(
            ElectrostaticPotentialSimulation_translation const&)
            = default;

    /**
     * @brief Get the electrostatic potential associated to the advection field.
     *
     * @param[in] coord
     *      The coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The electrostatic potential at the given coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordXY const coord, double const t) const
    {
        double const vx = ddc::get<X>(m_velocity);
        double const vy = ddc::get<Y>(m_velocity);
        return vy * ddc::get<X>(coord) - vx * ddc::get<Y>(coord);
    }

    /**
     * @brief Get the characteristic feet in the physical domain.
     *
     * @param[in] coord
     *      The original coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The characteristic feet in the physical domain.
     */
    KOKKOS_FUNCTION CoordXY exact_feet(CoordXY coord, double const t) const
    {
        return coord - t * m_velocity;
    }
};



/**
 * @brief Electrostatic potential for a rotation test of the 2D polar advection operator.
 *
 * The test advection field for a rotation in the physical domain is given by :
 *
 * @f$ A(x,y) = J_{\mathcal{F}}[v^r, v^\theta] = [-\partial_y \phi, \partial_x \phi]^T @f$,
 *
 * with @f$ v^r @f$ and @f$ v^\theta @f$ constants composants in the contravariant basis 
 * and @f$\mathcal{F}@f$ the circular mapping.
 * 
 * An adapted electrostatic potential for  @f$ v^r = 0 @f$ is 
 * @f$ \phi(x,y) = 1/2 v^\theta (x^2 + y^2 )@f$.
 *
 * The characteristic feet are then given by
 * - @f$ (X(t + dt), Y(t + dt))  = \mathcal{F} (R(t + dt), \Theta(t + dt))@f$,
 * - with
 *      - @f$ R(t + dt) = R(t) + dt v_r @f$,
 *      - @f$ \Theta(t + dt) = \Theta(t) + dt v_\theta @f$,
 *      - and @f$ (R(t), \Theta(t)) = \mathcal{F}^{-1} (X(t), Y(t))@f$.
 *
 */
class ElectrostaticPotentialSimulation_rotation
{
private:
    double const m_vr;
    double const m_vtheta;
    CartesianToCircular<X, Y, R, Theta> const m_mapping;

public:
    /**
     * @brief Instantiate an ElectrostaticPotentialSimulation_rotation advection field.
     *
     * @param[in] vtheta
     *      The constant second polar component of the advection field in the physical domain.
     */
    explicit ElectrostaticPotentialSimulation_rotation(CoordVtheta vtheta)
        : m_vr(0)
        , m_vtheta(vtheta)
        , m_mapping()
    {
    }

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION ElectrostaticPotentialSimulation_rotation(
            ElectrostaticPotentialSimulation_rotation const&)
            = default;

    /**
     * @brief Get the electrostatic potential associated to the advection field.
     *
     * @param[in] coord
     *      The coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The electrostatic potential at the given coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordXY const coord, double const t) const
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        return 0.5 * m_vtheta * (x * x + y * y);
    }

    /**
     * @brief Get the characteristic feet in the physical domain.
     *
     * @param[in] coord_xy
     *      The original coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The characteristic feet in the physical domain.
     */
    KOKKOS_FUNCTION CoordXY exact_feet(CoordXY coord_xy, double const t) const
    {
        CoordRTheta const coord_rtheta(m_mapping(coord_xy));
        CoordRTheta const velocity(m_vr, m_vtheta);
        CircularToCartesian<R, Theta, X, Y> logical_to_physical_mapping;
        return logical_to_physical_mapping(coord_rtheta - t * velocity);
    }
};



// TEST SIMULATIONS ------------------------------------------------------------------------------
/**
 * @brief Base class for the tests simulation of the 2D polar advection operator.
 *
 * The simulations are:
 * - TranslationAdvectionFieldSimulation,
 * - RotationAdvectionFieldSimulation,
 * - DecentredRotationAdvectionFieldSimulation.
 *
 * @see BslAdvectionPolar
 * @see FunctionToBeAdvected
 * @see ElectrostaticPotentialSimulation
 * @see AdvectionField
 */
template <class ElectrostaticPotentialSimulation, class FunctionToBeAdvected, class AdvectionField>
struct AdvectionFieldSimulation
{
    /// @brief The chosen electrostatical potential for the simulation.
    ElectrostaticPotentialSimulation const electrostatical_potential;

    ///@brief The chosen function to be advected for the simulation.
    FunctionToBeAdvected const function;

    ///@brief The chosen advection field for the simulation.
    AdvectionField const advection_field;
};


/**
 * @brief Simulation of a translated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a translation advection field defined by ElectrostaticPotentialSimulation_translation:
 *
 * - @f$ f(x, y) =  C \exp\left( - \frac{(x - x_0)^2}{2\sigma_x^2}  - \frac{(y - y_0)^2}{2\sigma_y^2}  \right)@f$
 * with
 *      - @f$ C = 1 @f$,
 *      - @f$ x_0 = -0.2 @f$ ,
 *      - @f$ y_0 = -0.2 @f$ ,
 *      - @f$ \sigma_x = 0.1 @f$ ,
 *      - @f$ \sigma_y = 0.1 @f$ ;
 *
 * and
 * - @f$ A(x,y) = [v_x, v_y]@f$,
 *
 * with @f$ v_x @f$ and @f$ v_y @f$ constants.
 *
 * @param[in] mapping
 *      The mapping from the logical domain to the physical domain.
 * @param[in] rmin
 *      The minimum value of @f$ r@f$ on the logical domain.
 * @param[in] rmax
 *      The maximum value of @f$ r@f$ on the logical domain.
 *
 * @see FunctionToBeAdvected_gaussian
 * @see ElectrostaticPotentialSimulation_translation
 * @see AdvectionField_translation
 */
template <class Mapping>
auto get_translation_advection_field_simulation(
        Mapping const& mapping,
        double const rmin,
        double const rmax)
{
    return AdvectionFieldSimulation<
            ElectrostaticPotentialSimulation_translation,
            FunctionToBeAdvected_gaussian<Mapping>,
            AdvectionField_translation<X, Y>>(
            {ElectrostaticPotentialSimulation_translation(
                     CoordVx(-std::cos(2 * M_PI * 511. / 4096.) / 2.),
                     CoordVy(-std::sin(2 * M_PI * 511. / 4096.) / 2.)),
             FunctionToBeAdvected_gaussian<Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax),
             AdvectionField_translation<X, Y>(DVector<X, Y>(
                     -std::cos(2 * M_PI * 511. / 4096.) / 2.,
                     -std::sin(2 * M_PI * 511. / 4096.) / 2.))});
}


/**
 * @brief Simulation of a rotated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a rotation advection field defined by ElectrostaticPotentialSimulation_rotation:
 *
 *
 * - @f$ f(x, y) =  C \exp\left( - \frac{(x - x_0)^2}{2\sigma_x^2}  - \frac{(y - y_0)^2}{2\sigma_y^2}  \right)@f$
 * with
 *      - @f$ C = 1 @f$,
 *      - @f$ x_0 = -0.2 @f$ ,
 *      - @f$ y_0 = -0.2 @f$ ,
 *      - @f$ \sigma_x = 0.1 @f$ ,
 *      - @f$ \sigma_y = 0.1 @f$ ;
 *
 * and
 * - @f$ A(x,y) = J_{\mathcal{F}}[v_r, v_\theta]@f$,
 *
 * with @f$ v_r @f$ and @f$ v_\theta @f$ constants and @f$ \mathcal{F}@f$ the
 * circular mapping.
 *
 * @param[in] mapping
 *      The mapping from the logical domain to the physical domain.
 * @param[in] rmin
 *      The minimum value of @f$ r@f$ on the logical domain.
 * @param[in] rmax
 *      The maximum value of @f$ r@f$ on the logical domain.
 *
 * @see FunctionToBeAdvected_gaussian
 * @see ElectrostaticPotentialSimulation_rotation
 * @see AdvectionField_rotation
 */
template <class Mapping>
auto get_rotation_advection_field_simulation(
        Mapping const& mapping,
        double const rmin,
        double const rmax)
{
    return AdvectionFieldSimulation<
            ElectrostaticPotentialSimulation_rotation,
            FunctionToBeAdvected_gaussian<Mapping>,
            AdvectionField_rotation<X, Y, R, Theta>>(
            {ElectrostaticPotentialSimulation_rotation(CoordVtheta(2 * M_PI)),
             FunctionToBeAdvected_gaussian<Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax),
             AdvectionField_rotation<X, Y, R, Theta>(DVector<R, Theta>(0, 2 * M_PI))});
}


/**
 * @brief Simulation of a decentred rotated ellipse-type function.
 *
 * The simulation is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * Define a simulation with a ellipse-type function defined by FunctionToBeAdvected_cos_4_ellipse
 * and a decentred rotation advection field defined by ElectrostaticPotentialSimulation_decentred_rotation:
 *
 * - @f$ f(x, y) =  \frac{1}{2}\left( G(r_1(x,y)) + G(r_2(x,y))\right)@f$
 *
 * with
 *      - @f$G(r) = \cos\left(\frac{\pi r}{2 a}\right)^4 * 1_{r<a}(r)@f$,
 *      - @f$r_1(x, y) = \sqrt{(x-x_0)^2 + 8(y-y_0)^2}@f$,
 *      - @f$r_2(x, y) = \sqrt{8(x-x_0)^2 + (y-y_0)^2}@f$,
 *
 * and
 * - @f$ A(x,y) = \omega [y_c - y, x - x_c]@f$,
 *
 * with @f$\omega = 2\pi@f$ and  @f$ (x_c, y_c) = (0.25, 0) @f$.
 *
 * @param[in] mapping
 *      The mapping from the logical domain to the physical domain.
 *
 * @see FunctionToBeAdvected_cos_4_ellipse
 * @see ElectrostaticPotentialSimulation_decentred_rotation
 * @see AdvectionField_decentred_rotation
 */
template <class Mapping>
auto get_decentred_rotation_advection_field_simulation(Mapping const& mapping)
{
    return AdvectionFieldSimulation<
            ElectrostaticPotentialSimulation_decentred_rotation,
            FunctionToBeAdvected_cos_4_ellipse<Mapping>,
            AdvectionField_decentred_rotation<X, Y>>(
            {ElectrostaticPotentialSimulation_decentred_rotation(),
             FunctionToBeAdvected_cos_4_ellipse<Mapping>(mapping),
             AdvectionField_decentred_rotation<X, Y>()});
}
