// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "../advection_2d_rtheta/test_cases.hpp"

#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "math_tools.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "spline_interpolator_2d_rp.hpp"

/*
 *  This file defines
 *   - the TEST ELECTROSTATICAL POTENTIALS,
 *   - and the TEST SIMULATIONS.
 */



// TEST ELECTROSTATICAL POTENTIALS ------------------------------------------------
/**
 * @brief Base class for the advection field for the tests of the 2D polar advection
 * operator.
 *
 * @see BslAdvectionRTheta
 * @see AdvectionSimulation
 */
class ElectrostaticalPotentialSimulation
{
public:
    virtual ~ElectrostaticalPotentialSimulation() = default;

    /**
     * @brief Get the advection field in the physical index range.
     *
     * @param[in] coord
     *      The coordinate in the physical index range.
     * @param[in] t
     *      Time component.
     *
     * @return The advection field in the physical index range.
     */
    virtual double operator()(CoordXY const coord, double const t = 0.) const = 0;

    /**
     * @brief Get the characteristic feet in the physical index range.
     *
     * @param[in] coord
     *      The original coordinate in the physical index range.
     * @param[in] t
     *      Time component.
     *
     * @return The characteristic feet in the physical index range.
     */
    virtual CoordXY exact_feet(CoordXY coord, double t) const = 0;
};


/**
 * @brief Advection field for the tests of the 2D polar advection operator.
 *
 *
 * The test advection field is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889):
 *
 * @f$ A(x,y) = \omega [y_c - y, x - x_c]@f$.
 *
 * The characteristic feet are then given by
 * - @f$ X(t + dt) = x_c + (X(t) - x_c) \cos(\omega dt) - (Y(t) - y_c) \sin(\omega dt) @f$,
 * - @f$ Y(t + dt) = y_c + (X(t) - x_c) \sin(\omega dt) + (Y(t) - y_c) \cos(\omega dt) @f$.
 *
 */
class ElectrostaticalPotentialSimulation_decentred_rotation
    : public ElectrostaticalPotentialSimulation
{
private:
    double const m_omega;
    double const m_xc;
    double const m_yc;
    double const m_x_bar;
    double const m_y_bar;

public:
    /**
     * @brief Instantiate an ElectrostaticalPotentialSimulation_decentred_rotation advection field.
     *
     */
    ElectrostaticalPotentialSimulation_decentred_rotation()
        : m_omega(2 * M_PI)
        , m_xc(0.25)
        , m_yc(0.)
        , m_x_bar(0.)
        , m_y_bar(0.)
    {
    }

    double operator()(CoordXY const coord, double const t) const final
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        return m_omega * (-0.5 * x * x + m_xc * x - 0.5 * y * y + m_yc * y);
    }

    CoordXY exact_feet(CoordXY coord, double const t) const final
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        double const foot_x
                = m_xc + (x - m_xc) * std::cos(m_omega * -t) - (y - m_yc) * std::sin(m_omega * -t);
        double const foot_y
                = m_yc + (x - m_xc) * std::sin(m_omega * -t) + (y - m_yc) * std::cos(m_omega * -t);
        return CoordXY(foot_x, foot_y);
    }
};



/**
 * @brief Advection field for the tests of the 2D polar advection operator.
 *
 *
 * The test advection field for a translation in the physical index range is given by :
 *
 * @f$ A(x,y) = [v_x, v_y]@f$.
 *
 * with @f$ v_x @f$ and @f$ v_y @f$ constants.
 *
 * The characteristic feet are then given by
 * - @f$ X(t + dt) = X(t) + dt v_x @f$,
 * - @f$ Y(t + dt) = Y(t) + dt v_y @f$.
 *
 */
class ElectrostaticalPotentialSimulation_translation : public ElectrostaticalPotentialSimulation
{
private:
    CoordXY const m_velocity;

public:
    /**
     * @brief Instantiate an ElectrostaticalPotentialSimulation_translation advection field.
     *
     * @param[in] vx
     *      The constant first component of the advection field in the physical index range.
     * @param[in] vy
     *      The constant second component of the advection field in the physical index range.
     */
    ElectrostaticalPotentialSimulation_translation(CoordVx vx, CoordVy vy)
        : m_velocity(ddc::get<Vx>(vx), ddc::get<Vy>(vy))
    {
    }

    double operator()(CoordXY const coord, double const t) const final
    {
        double const vx = ddc::get<X>(m_velocity);
        double const vy = ddc::get<Y>(m_velocity);
        return -vy * ddc::get<X>(coord) + vx * ddc::get<Y>(coord);
    }

    CoordXY exact_feet(CoordXY coord, double const t) const final
    {
        return coord - t * m_velocity;
    }
};



/**
 * @brief Advection field for the tests of the 2D polar advection operator.
 *
 *
 * The test advection field for a rotation in the physical index range is given by :
 *
 * @f$ A(x,y) = J_{\mathcal{F}}[v_r, v_\theta]@f$.
 *
 * with @f$ v_r @f$ and @f$ v_\theta @f$ constants and @f$ \mathcal{F}@f$ the
 * circular mapping.
 *
 * The characteristic feet are then given by
 * - @f$ (X(t + dt), Y(t + dt))  = \mathcal{F} (R(t + dt), \Theta(t + dt))@f$,
 * - with
 *      - @f$ R(t + dt) = R(t) + dt v_r @f$,
 *      - @f$ \Theta(t + dt) = \Theta(t) + dt v_\theta @f$,
 *      - and @f$ (R(t), \Theta(t)) = \mathcal{F}^{-1} (X(t), Y(t))@f$.
 *
 */
class ElectrostaticalPotentialSimulation_rotation : public ElectrostaticalPotentialSimulation
{
private:
    double const m_vr;
    double const m_vtheta;
    CartesianToCircular<X, Y, R, Theta> const m_mapping;

public:
    /**
     * @brief Instantiate an ElectrostaticalPotentialSimulation_rotation advection field.
     *
     * @param[in] vtheta
     *      The constant second polar component of the advection field in the physical index range.
     */
    explicit ElectrostaticalPotentialSimulation_rotation(CoordVtheta vtheta)
        : m_vr(0)
        , m_vtheta(vtheta)
        , m_mapping()
    {
    }

    double operator()(CoordXY const coord, double const t) const final
    {
        CoordRTheta const coord_rp(m_mapping(coord));
        double const r = ddc::get<R>(coord_rp);
        return -0.5 * r * r * m_vtheta;
    }

    CoordXY exact_feet(CoordXY coord_xy, double const t) const final
    {
        CoordRTheta const coord_rp(m_mapping(coord_xy));
        CoordRTheta const velocity(m_vr, m_vtheta);
        CircularToCartesian<R, Theta, X, Y> logical_to_physical_mapping;
        return logical_to_physical_mapping(coord_rp - t * velocity);
    }
};



// TEST SIMULATIONS ---------------------------------------------------------------
/**
 * @brief Base class for the tests simulation of the 2D polar advection operator.
 *
 * The simulations are:
 * - TranslationAdvectionFieldSimulation,
 * - RotationAdvectionFieldSimulation,
 * - DecentredRotationAdvectionFieldSimulation.
 *
 * @see BslAdvectionRTheta
 * @see FunctionToBeAdvected
 * @see ElectrostaticalPotentialSimulation
 * @see AdvectionField
 */
template <
        class ElectrostaticalPotentialSimulation,
        class FunctionToBeAdvected,
        class AdvectionField>
class AdvectionFieldSimulation
{
protected:
    /**
     * @brief The chosen electrostatical potential for the simulation.
     */
    ElectrostaticalPotentialSimulation const m_electrostatical_potential;
    /**
     * @brief The chosen function to be advected for the simulation.
     */
    FunctionToBeAdvected const m_function;
    /**
     * @brief The chosen advection field for the simulation.
     */
    AdvectionField const m_advection_field;


public:
    /**
     * @brief Instantiate a AdvectionSimulation simulation.
     *
     * @param[in] electrostatical_potential
     *      An ElectrostaticalPotentialSimulation type object.
     * @param[in] function
     *      A FunctionToBeAdvected type object.
     * @param[in] advection_field
     *      A AdvectionField type object.
     */
    AdvectionFieldSimulation(
            ElectrostaticalPotentialSimulation const electrostatical_potential,
            FunctionToBeAdvected const function,
            AdvectionField const advection_field)
        : m_electrostatical_potential(electrostatical_potential)
        , m_function(function)
        , m_advection_field(advection_field)
    {
    }
    virtual ~AdvectionFieldSimulation() = default;


    /**
     * @brief Get the electrostatical potential of the simulation.
     *
     * @return A constant reference to the electrostatical potential created in 
     * the ElectrostaticalPotentialSimulation child class.
     */
    ElectrostaticalPotentialSimulation const& get_electrostatique_potential() const
    {
        return m_electrostatical_potential;
    }

    /**
     * @brief Get the test function of the simulation.
     *
     * @return A constant reference to the test function created in the AdvectionSimulation child class.
     */
    FunctionToBeAdvected const& get_test_function() const
    {
        return m_function;
    }


    /**
     * @brief Get the advection field of the simulation.
     *
     * @return  A constant reference to the advection field created in the AdvectionSimulation child class.
     */
    AdvectionField const& get_advection_field() const
    {
        return m_advection_field;
    }
};


/**
 * @brief Simulation of a translated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a translation advection field defined by ElectrostaticalPotentialSimulation_translation:
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
 * @see FunctionToBeAdvected_gaussian
 * @see ElectrostaticalPotentialSimulation_translation
 * @see AdvectionField_translation
 */
template <class Mapping>
class TranslationAdvectionFieldSimulation
    : public AdvectionFieldSimulation<
              ElectrostaticalPotentialSimulation_translation,
              FunctionToBeAdvected_gaussian<Mapping>,
              AdvectionField_translation>
{
public:
    /**
     * @brief Instantiate a TranslationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range.
     * @param[in] rmin
     *      The minimum value of @f$ r@f$ on the logical index range.
     * @param[in] rmax
     *      The maximum value of @f$ r@f$ on the logical index range.
     */
    TranslationAdvectionFieldSimulation(
            Mapping const& mapping,
            double const rmin,
            double const rmax)
        : AdvectionFieldSimulation<
                ElectrostaticalPotentialSimulation_translation,
                FunctionToBeAdvected_gaussian<Mapping>,
                AdvectionField_translation>(
                ElectrostaticalPotentialSimulation_translation(
                        CoordVx(-std::cos(2 * M_PI * 511. / 4096.) / 2.),
                        CoordVy(-std::sin(2 * M_PI * 511. / 4096.) / 2.)),
                FunctionToBeAdvected_gaussian<
                        Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax),
                AdvectionField_translation(
                        CoordVx(-std::cos(2 * M_PI * 511. / 4096.) / 2.),
                        CoordVy(-std::sin(2 * M_PI * 511. / 4096.) / 2.)))
    {
    }
};


/**
 * @brief Simulation of a rotated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a rotation advection field defined by ElectrostaticalPotentialSimulation_rotation:
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
 * @see FunctionToBeAdvected_gaussian
 * @see ElectrostaticalPotentialSimulation_rotation
 * @see AdvectionField_rotation
 */
template <class Mapping>
class RotationAdvectionFieldSimulation
    : public AdvectionFieldSimulation<
              ElectrostaticalPotentialSimulation_rotation,
              FunctionToBeAdvected_gaussian<Mapping>,
              AdvectionField_rotation>
{
public:
    /**
     * @brief Instantiate a RotationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range.
     * @param[in] rmin
     *      The minimum value of @f$ r@f$ on the logical index range.
     * @param[in] rmax
     *      The maximum value of @f$ r@f$ on the logical index range.
     */
    RotationAdvectionFieldSimulation(Mapping const& mapping, double const rmin, double const rmax)
        : AdvectionFieldSimulation<
                ElectrostaticalPotentialSimulation_rotation,
                FunctionToBeAdvected_gaussian<Mapping>,
                AdvectionField_rotation>(
                ElectrostaticalPotentialSimulation_rotation(CoordVtheta(2 * M_PI)),
                FunctionToBeAdvected_gaussian<
                        Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax),
                AdvectionField_rotation(CoordVr(0), CoordVtheta(2 * M_PI)))
    {
    }
};


/**
 * @brief Simulation of a decentred rotated elipse-type function.
 *
 * The simulation is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * Define a simulation with a elipse-type function defined by FunctionToBeAdvected_cos_4_elipse
 * and a decentred rotation advection field defined by ElectrostaticalPotentialSimulation_decentred_rotation:
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
 *
 * @see FunctionToBeAdvected_cos_4_elipse
 * @see ElectrostaticalPotentialSimulation_decentred_rotation
 * @see AdvectionField_decentred_rotation
 */
template <class Mapping>
class DecentredRotationAdvectionFieldSimulation
    : public AdvectionFieldSimulation<
              ElectrostaticalPotentialSimulation_decentred_rotation,
              FunctionToBeAdvected_cos_4_elipse<Mapping>,
              AdvectionField_decentred_rotation>
{
public:
    /**
     * @brief Instantiate a DecentredRotationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range.
     */
    explicit DecentredRotationAdvectionFieldSimulation(Mapping const& mapping)
        : AdvectionFieldSimulation<
                ElectrostaticalPotentialSimulation_decentred_rotation,
                FunctionToBeAdvected_cos_4_elipse<Mapping>,
                AdvectionField_decentred_rotation>(
                ElectrostaticalPotentialSimulation_decentred_rotation(),
                FunctionToBeAdvected_cos_4_elipse<Mapping>(mapping),
                AdvectionField_decentred_rotation())
    {
    }
};
