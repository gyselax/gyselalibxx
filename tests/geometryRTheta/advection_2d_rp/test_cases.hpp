// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "math_tools.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "spline_interpolator_2d_rp.hpp"

/*
 *  This file defines
 *   - the TEST FUNCTIONS,
 *   - the TEST ADVECTION FIELDS,
 *   - and the TEST SIMULATIONS.
 */


//---------------------------------------------------------------------------------
//                           TEST FUNCTIONS
//---------------------------------------------------------------------------------

/**
 * @brief Test function for the 2D polar advection operator.
 *
 * The test function is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889):
 *
 * @f$ f(x, y) =  \frac{1}{2}\left( G(r_1(x,y)) + G(r_2(x,y))\right)@f$
 *
 * with
 *  - @f$G(r) = \cos\left(\frac{\pi r}{2 a}\right)^4 * 1_{r<a}(r)@f$,
 *  - @f$r_1(x, y) = \sqrt{(x-x_0)^2 + 8(y-y_0)^2}@f$,
 *  - @f$r_2(x, y) = \sqrt{8(x-x_0)^2 + (y-y_0)^2}@f$.
 *
 */
template <class Mapping>
class FunctionToBeAdvected_cos_4_ellipse
{
private:
    Mapping const& m_mapping;

public:
    /**
     * @brief Instantiate a FunctionToBeAdvected_cos_4_ellipse function.
     *
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range.
     */
    explicit FunctionToBeAdvected_cos_4_ellipse(Mapping const& mapping) : m_mapping(mapping) {}

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION FunctionToBeAdvected_cos_4_ellipse(
            FunctionToBeAdvected_cos_4_ellipse const&)
            = default;

    /**
     * @brief Get the value of the function.
     *
     * @param[in] coord_rp
     *      The coordinate where we want to evaluate the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta coord_rp) const
    {
        CoordXY const coord_xy(m_mapping(coord_rp));
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);

        double const a = 0.3;
        double const x_bar = 0.;
        double const y_bar = 0.;

        double const r_1 = Kokkos::sqrt((x - x_bar) * (x - x_bar) + 8 * (y - y_bar) * (y - y_bar));
        double const r_2 = Kokkos::sqrt(8 * (x - x_bar) * (x - x_bar) + (y - y_bar) * (y - y_bar));

        double const G_1 = ipow(Kokkos::cos(M_PI * r_1 / 2. / a), 4) * (r_1 < a);
        double const G_2 = ipow(Kokkos::cos(M_PI * r_2 / 2. / a), 4) * (r_2 < a);

        return 1. / 2. * (G_1 + G_2);
    }
};


/**
 * @brief Test function for the 2D polar advection operator.
 *
 * The test function is a Gaussian centred in @f$ (x_0, y_0 @f$):
 *
 * @f$ f(x, y) =  C \exp\left( - \frac{(x - x_0)^2}{2\sigma_x^2}  - \frac{(y - y_0)^2}{2\sigma_y^2}  \right)@f$
 *
 */
template <class Mapping>
class FunctionToBeAdvected_gaussian
{
private:
    Mapping const& m_mapping;
    double const m_constant;
    double const m_x0;
    double const m_y0;
    double const m_sig_x;
    double const m_sig_y;
    double const m_rmin;
    double const m_rmax;

public:
    /**
     * @brief Instantiate a FunctionToBeAdvected_gaussian function.
     *
     * @param[in] mapping
     *      The mapping from the logical index range to the physical index range.
     * @param[in] constant
     *      The constant @f$ C@f$ in front of the exponential.
     * @param[in] x0
     *      The first component @f$x_0@f$ of the point where the Gaussian is centred.
     * @param[in] y0
     *      The second component @f$y_0@f$ of the point where the Gaussian is centred.
     * @param[in] sig_x
     *      The standard deviation @f$\sigma_x @f$ in the first dimension.
     * @param[in] sig_y
     *      The standard deviation @f$\sigma_y @f$ in the second dimension.
     * @param[in] rmin
     *      The minimum value of @f$ r@f$.
     * @param[in] rmax
     *      The maximum value of @f$ r@f$.
     */
    FunctionToBeAdvected_gaussian(
            Mapping const& mapping,
            double const constant,
            double const x0,
            double const y0,
            double const sig_x,
            double const sig_y,
            double const rmin,
            double const rmax)
        : m_mapping(mapping)
        , m_constant(constant)
        , m_x0(x0)
        , m_y0(y0)
        , m_sig_x(sig_x)
        , m_sig_y(sig_y)
        , m_rmin(rmin)
        , m_rmax(rmax)
    {
    }

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION FunctionToBeAdvected_gaussian(
            FunctionToBeAdvected_gaussian const&)
            = default;

    /**
     * @brief Get the value of the function.
     *
     * @param[in] coord_rp
     *      The coordinate where we want to evaluate the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordRTheta coord_rp) const
    {
        // Gaussian centred in (x0, y0):
        CoordXY const coord_xy(m_mapping(coord_rp));
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);
        double const r = ddc::get<R>(coord_rp);
        if ((m_rmin <= r) and (r <= m_rmax)) {
            return m_constant
                   * Kokkos::exp(
                           -ipow(x - m_x0, 2) / (2 * m_sig_x * m_sig_x)
                           - ipow(y - m_y0, 2) / (2 * m_sig_y * m_sig_y));
        } else {
            return 0.0;
        }
    }
};



//---------------------------------------------------------------------------------
//                           TEST ADVECTION FIELDS
//---------------------------------------------------------------------------------

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
class AdvectionField_decentred_rotation
{
private:
    double const m_omega;
    double const m_xc;
    double const m_yc;
    double const m_x_bar;
    double const m_y_bar;

public:
    /**
     * @brief Instantiate an AdvectionField_decentred_rotation advection field.
     *
     */
    AdvectionField_decentred_rotation()
        : m_omega(2 * M_PI)
        , m_xc(0.25)
        , m_yc(0.)
        , m_x_bar(0.)
        , m_y_bar(0.)
    {
    }

    /// Copy operator
    explicit KOKKOS_DEFAULTED_FUNCTION AdvectionField_decentred_rotation(
            AdvectionField_decentred_rotation const&)
            = default;

    /**
￼     * @brief Get the advection field in the physical index range.
￼     *
￼     * @param[in] coord
￼     *      The coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The advection field in the physical index range.
￼     */
    KOKKOS_FUNCTION CoordXY operator()(CoordXY const coord, double const t) const
    {
        double const x = m_omega * (m_yc - ddc::get<Y>(coord));
        double const y = m_omega * (ddc::get<X>(coord) - m_xc);
        return CoordXY(x, y);
    }

    /**
￼     * @brief Get the characteristic feet in the physical index range.
￼     *
￼     * @param[in] coord
￼     *      The original coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The characteristic feet in the physical index range.
￼     */
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
class AdvectionField_translation
{
private:
    CoordXY const m_velocity;

public:
    /**
     * @brief Instantiate an AdvectionField_translation advection field.
     *
     * @param[in] vx
     *      The constant first component of the advection field in the physical index range.
     * @param[in] vy
     *      The constant second component of the advection field in the physical index range.
     */
    AdvectionField_translation(CoordVx vx, CoordVy vy)
        : m_velocity(ddc::get<Vx>(vx), ddc::get<Vy>(vy))
    {
    }

    /// Copy operator
    KOKKOS_DEFAULTED_FUNCTION AdvectionField_translation(AdvectionField_translation const&)
            = default;

    /**
￼     * @brief Get the advection field in the physical index range.
￼     *
￼     * @param[in] coord
￼     *      The coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The advection field in the physical index range.
￼     */
    KOKKOS_FUNCTION CoordXY operator()(CoordXY const coord, double const t) const
    {
        return m_velocity;
    }

    /**
￼     * @brief Get the characteristic feet in the physical index range.
￼     *
￼     * @param[in] coord
￼     *      The original coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The characteristic feet in the physical index range.
￼     */
    KOKKOS_FUNCTION CoordXY exact_feet(CoordXY coord, double const t) const
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
class AdvectionField_rotation
{
private:
    double const m_vr;
    double const m_vtheta;
    CartesianToCircular<X, Y, R, Theta> const m_physical_to_logical_mapping;
    CircularToCartesian<R, Theta, X, Y> const m_logical_to_physical_mapping;


public:
    /**
     * @brief Instantiate an AdvectionField_rotation advection field.
     *
     * @param[in] vr
     *      The constant first polar component of the advection field in the physical index range.
     * @param[in] vtheta
     *      The constant second polar component of the advection field in the physical index range.
     */
    AdvectionField_rotation(CoordVr vr, CoordVtheta vtheta)
        : m_vr(vr)
        , m_vtheta(vtheta)
        , m_physical_to_logical_mapping()
        , m_logical_to_physical_mapping()
    {
    }

    /// Copy operator
    KOKKOS_DEFAULTED_FUNCTION AdvectionField_rotation(AdvectionField_rotation const&) = default;

    /**
￼     * @brief Get the advection field in the physical index range.
￼     *
￼     * @param[in] coord
￼     *      The coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The advection field in the physical index range.
￼     */
    KOKKOS_FUNCTION CoordXY operator()(CoordXY const coord, double const t) const
    {
        CoordRTheta const coord_rp(m_physical_to_logical_mapping(coord));
        std::array<std::array<double, 2>, 2> jacobian;
        m_logical_to_physical_mapping.jacobian_matrix(coord_rp, jacobian);
        double const vx = m_vr * jacobian[0][0] + m_vtheta * jacobian[0][1];
        double const vy = m_vr * jacobian[1][0] + m_vtheta * jacobian[1][1];
        return CoordXY(vx, vy);
    }

    /**
￼     * @brief Get the characteristic feet in the physical index range.
￼     *
￼     * @param[in] coord_xy
￼     *      The original coordinate in the physical index range.
￼     * @param[in] t
￼     *      Time component.
￼     *
￼     * @return The characteristic feet in the physical index range.
￼     */
    KOKKOS_FUNCTION CoordXY exact_feet(CoordXY coord_xy, double const t) const
    {
        CoordRTheta const coord_rp(m_physical_to_logical_mapping(coord_xy));
        CoordRTheta const velocity(m_vr, m_vtheta);
        return m_logical_to_physical_mapping(coord_rp - t * velocity);
    }
};



// TEST SIMULATIONS ---------------------------------------------------------------
/**
 * @brief Base class for the tests simulation of the 2D polar advection operator.
 *
 * The simulations are:
 * - TranslationSimulation,
 * - RotationSimulation,
 * - DecentredRotationSimulation.
 *
 * @see BslAdvectionRTheta
 * @see FunctionToBeAdvected
 * @see AdvectionField
 */
template <class AdvectionField, class FunctionToBeAdvected>
struct AdvectionSimulation
{
    /**
     * @brief The chosen advection field for the simulation.
     */
    AdvectionField const advection_field;
    /**
     * @brief The chosen function to be advected for the simulation.
     */
    FunctionToBeAdvected const advected_function;
};

/**
 * @brief Simulation of a translated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a translation advection field defined by AdvectionField_translation:
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
 * @see AdvectionField_translation
 */
template <class Mapping>
AdvectionSimulation<AdvectionField_translation, FunctionToBeAdvected_gaussian<Mapping>>
get_translation_simulation(Mapping const& mapping, double const rmin, double const rmax)
{
    return AdvectionSimulation<AdvectionField_translation, FunctionToBeAdvected_gaussian<Mapping>>(
            {AdvectionField_translation(
                     CoordVx(std::cos(2 * M_PI * 511. / 4096.) / 2.),
                     CoordVy(std::sin(2 * M_PI * 511. / 4096.) / 2.)),
             FunctionToBeAdvected_gaussian<
                     Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax)});
}

/**
 * @brief Simulation of a rotated Gaussian.
 *
 * Define a simulation with a Gaussian defined by FunctionToBeAdvected_gaussian
 * and a rotation advection field defined by AdvectionField_rotation:
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
 * @see AdvectionField_rotation
 */
template <class Mapping>
AdvectionSimulation<AdvectionField_rotation, FunctionToBeAdvected_gaussian<Mapping>>
get_rotation_simulation(Mapping const& mapping, double const rmin, double const rmax)
{
    return AdvectionSimulation<AdvectionField_rotation, FunctionToBeAdvected_gaussian<Mapping>>(
            {AdvectionField_rotation(CoordVr(0.), CoordVtheta(2 * M_PI)),
             FunctionToBeAdvected_gaussian<
                     Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax)});
}

/**
 * @brief Simulation of a decentred rotated ellipse-type function.
 *
 * The simulation is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * Define a simulation with a ellipse-type function defined by FunctionToBeAdvected_cos_4_ellipse
 * and a decentred rotation advection field defined by AdvectionField_decentred_rotation:
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
 * @see FunctionToBeAdvected_cos_4_ellipse
 * @see AdvectionField_decentred_rotation
 */
template <class Mapping>
AdvectionSimulation<AdvectionField_decentred_rotation, FunctionToBeAdvected_cos_4_ellipse<Mapping>>
get_decentred_rotation_simulation(Mapping const& mapping)
{
    return AdvectionSimulation<
            AdvectionField_decentred_rotation,
            FunctionToBeAdvected_cos_4_ellipse<Mapping>>(
            {AdvectionField_decentred_rotation(),
             FunctionToBeAdvected_cos_4_ellipse<Mapping>(mapping)});
}
