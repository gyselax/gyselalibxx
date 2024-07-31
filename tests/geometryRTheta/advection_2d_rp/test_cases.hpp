#pragma once

#include <ddc/ddc.hpp>

#include <sll/gauss_legendre_integration.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/math_tools.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "spline_interpolator_2d_rp.hpp"

/*
 *  This file defines
 *   - the TEST FUNCTIONS,
 *   - the TEST ADVECTION FIELDS,
 *   - and the TEST SIMULATIONS.
 */



// TEST FUNCTIONS -----------------------------------------------------------------
/**
 * @brief Base class for the test functions for the 2D polar advection operator.
 *
 * @see BslAdvectionRP
 * @see AdvectionSimulation
 */
class FunctionToBeAdvected
{
public:
    virtual ~FunctionToBeAdvected() = default;

    /**
     * @brief Get the value of the function.
     *
     * @param[in] coord
     *      The coordinate where we want to evaluate the function.
     *
     * @return The value of the function at the coordinate.
     */
    virtual double operator()(CoordRP coord) = 0;
};


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
class FunctionToBeAdvected_cos_4_elipse : public FunctionToBeAdvected
{
private:
    Mapping const& m_mapping;

public:
    /**
     * @brief Instantiate a FunctionToBeAdvected_cos_4_elipse function.
     *
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain.
     */
    FunctionToBeAdvected_cos_4_elipse(Mapping const& mapping) : m_mapping(mapping) {};
    ~FunctionToBeAdvected_cos_4_elipse() {};

    double operator()(CoordRP coord_rp) final
    {
        CoordXY const coord_xy(m_mapping(coord_rp));
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);

        double const a = 0.3;
        double const x_bar = 0.;
        double const y_bar = 0.;

        double const r_1 = std::sqrt((x - x_bar) * (x - x_bar) + 8 * (y - y_bar) * (y - y_bar));
        double const r_2 = std::sqrt(8 * (x - x_bar) * (x - x_bar) + (y - y_bar) * (y - y_bar));

        double const G_1 = std::pow(std::cos(M_PI * r_1 / 2. / a), 4) * (r_1 < a);
        double const G_2 = std::pow(std::cos(M_PI * r_2 / 2. / a), 4) * (r_2 < a);

        return 1. / 2. * (G_1 + G_2);
    };
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
class FunctionToBeAdvected_gaussian : public FunctionToBeAdvected
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
     *      The mapping from the logical domain to the physical domain.
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
        , m_rmax(rmax) {};
    ~FunctionToBeAdvected_gaussian() {};

    double operator()(CoordRP coord_rp) final
    {
        // Gaussian centered in (x0, y0):
        CoordXY const coord_xy(m_mapping(coord_rp));
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);
        double const r = ddc::get<RDimR>(coord_rp);
        if ((m_rmin <= r) and (r <= m_rmax)) {
            return m_constant
                   * std::exp(
                           -pow(x - m_x0, 2) / (2 * m_sig_x * m_sig_x)
                           - pow(y - m_y0, 2) / (2 * m_sig_y * m_sig_y));
        } else {
            return 0.0;
        }
    };
};



// TEST ADVECTION FIELDS ----------------------------------------------------------
/**
 * @brief Base class for the advection field for the tests of the 2D polar advection
 * operator.
 *
 * @see BslAdvectionRP
 * @see AdvectionSimulation
 */
class AdvectionField
{
public:
    virtual ~AdvectionField() = default;

    /**
     * @brief Get the advection field in the physical domain.
     *
     * @param[in] coord
     *      The coordinate in the physical domain.
     * @param[in] t
     *      Time component.
     *
     * @return The advection field in the physical domain.
     */
    virtual CoordXY operator()(CoordXY const coord, double const t = 0.) const = 0;

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
class AdvectionField_decentred_rotation : public AdvectionField
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
        , m_y_bar(0.) {};
    ~AdvectionField_decentred_rotation() {};

    CoordXY operator()(CoordXY const coord, double const t) const final
    {
        double const x = m_omega * (m_yc - ddc::get<RDimY>(coord));
        double const y = m_omega * (ddc::get<RDimX>(coord) - m_xc);
        return CoordXY(x, y);
    };

    CoordXY exact_feet(CoordXY coord, double const t) const final
    {
        double const x = ddc::get<RDimX>(coord);
        double const y = ddc::get<RDimY>(coord);
        double const foot_x
                = m_xc + (x - m_xc) * std::cos(m_omega * -t) - (y - m_yc) * std::sin(m_omega * -t);
        double const foot_y
                = m_yc + (x - m_xc) * std::sin(m_omega * -t) + (y - m_yc) * std::cos(m_omega * -t);
        return CoordXY(foot_x, foot_y);
    };
};



/**
 * @brief Advection field for the tests of the 2D polar advection operator.
 *
 *
 * The test advection field for a translation in the physical domain is given by :
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
class AdvectionField_translation : public AdvectionField
{
private:
    CoordXY const m_velocity;

public:
    /**
     * @brief Instantiate an AdvectionField_translation advection field.
     *
     * @param[in] vx
     *      The constant first component of the advection field in the physical domain.
     * @param[in] vy
     *      The constant second component of the advection field in the physical domain.
     */
    AdvectionField_translation(CoordVx vx, CoordVy vy)
        : m_velocity(ddc::get<RDimVx>(vx), ddc::get<RDimVy>(vy)) {};
    ~AdvectionField_translation() {};

    CoordXY operator()(CoordXY const coord, double const t) const final
    {
        return m_velocity;
    };

    CoordXY exact_feet(CoordXY coord, double const t) const final
    {
        return coord - t * m_velocity;
    };
};



/**
 * @brief Advection field for the tests of the 2D polar advection operator.
 *
 *
 * The test advection field for a rotation in the physical domain is given by :
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
class AdvectionField_rotation : public AdvectionField
{
private:
    double const m_vr;
    double const m_vp;
    CircularToCartesian<RDimX, RDimY, RDimR, RDimP> const m_mapping;

public:
    /**
     * @brief Instantiate an AdvectionField_rotation advection field.
     *
     * @param[in] vr
     *      The constant first polar component of the advection field in the physical domain.
     * @param[in] vp
     *      The constant second polar component of the advection field in the physical domain.
     */
    AdvectionField_rotation(CoordVr vr, CoordVp vp) : m_vr(vr), m_vp(vp), m_mapping() {};
    ~AdvectionField_rotation() {};

    CoordXY operator()(CoordXY const coord, double const t) const final
    {
        CoordRP const coord_rp(m_mapping(coord));
        std::array<std::array<double, 2>, 2> jacobian;
        m_mapping.jacobian_matrix(coord_rp, jacobian);
        double const vx = m_vr * jacobian[0][0] + m_vp * jacobian[0][1];
        double const vy = m_vr * jacobian[1][0] + m_vp * jacobian[1][1];
        return CoordXY(vx, vy);
    };

    CoordXY exact_feet(CoordXY coord_xy, double const t) const final
    {
        CoordRP const coord_rp(m_mapping(coord_xy));
        CoordRP const velocity(m_vr, m_vp);
        return m_mapping(coord_rp - t * velocity);
    };
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
 * @see BslAdvectionRP
 * @see FunctionToBeAdvected
 * @see AdvectionField
 */
template <class AdvectionField, class FunctionToBeAdvected>
class AdvectionSimulation
{
protected:
    /**
     * @brief The chosen advection field for the simulation.
     */
    AdvectionField const m_advection_field;
    /**
     * @brief The chosen function to be advected for the simulation.
     */
    FunctionToBeAdvected const m_function;

public:
    /**
     * @brief Instantiate a AdvectionSimulation simulation.
     *
     * @param[in] advection_field
     *      An AdvectionField type object.
     * @param[in] function
     *      A FunctionToBeAdvected type object.
     */
    AdvectionSimulation(AdvectionField const advection_field, FunctionToBeAdvected const function)
        : m_advection_field(advection_field)
        , m_function(function) {};
    virtual ~AdvectionSimulation() = default;


    /**
     * @brief Get the advection field of the simulation.
     *
     * @return A constant reference to the advection field created in the AdvectionSimulation child class.
     */
    AdvectionField const& get_advection_field() const
    {
        return m_advection_field;
    };

    /**
     * @brief Get the test function of the simulation.
     *
     * @return A constant reference to the test function created in the AdvectionSimulation child class.
     */
    FunctionToBeAdvected const& get_test_function() const
    {
        return m_function;
    };
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
class TranslationSimulation
    : public AdvectionSimulation<AdvectionField_translation, FunctionToBeAdvected_gaussian<Mapping>>
{
public:
    /**
     * @brief Instantiate a TranslationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] rmin
     *      The minimum value of @f$ r@f$ on the logical domain.
     * @param[in] rmax
     *      The maximum value of @f$ r@f$ on the logical domain.
     */
    TranslationSimulation(Mapping const& mapping, double const rmin, double const rmax)
        : AdvectionSimulation<AdvectionField_translation, FunctionToBeAdvected_gaussian<Mapping>>(
                AdvectionField_translation(
                        CoordVx(std::cos(2 * M_PI * 511. / 4096.) / 2.),
                        CoordVy(std::sin(2 * M_PI * 511. / 4096.) / 2.)),
                FunctionToBeAdvected_gaussian<
                        Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax)) {};
};


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
class RotationSimulation
    : public AdvectionSimulation<AdvectionField_rotation, FunctionToBeAdvected_gaussian<Mapping>>
{
public:
    /**
     * @brief Instantiate a RotationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain.
     * @param[in] rmin
     *      The minimum value of @f$ r@f$ on the logical domain.
     * @param[in] rmax
     *      The maximum value of @f$ r@f$ on the logical domain.
     */
    RotationSimulation(Mapping const& mapping, double const rmin, double const rmax)
        : AdvectionSimulation<AdvectionField_rotation, FunctionToBeAdvected_gaussian<Mapping>>(
                AdvectionField_rotation(CoordVr(0.), CoordVp(2 * M_PI)),
                FunctionToBeAdvected_gaussian<
                        Mapping>(mapping, 1., -0.2, -0.2, 0.1, 0.1, rmin, rmax)) {};
};


/**
 * @brief Simulation of a decentred rotated elipse-type function.
 *
 * The simulation is given in Edoardo Zoni's article
 * (https://doi.org/10.1016/j.jcp.2019.108889).
 *
 * Define a simulation with a elipse-type function defined by FunctionToBeAdvected_cos_4_elipse
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
 * @see FunctionToBeAdvected_cos_4_elipse
 * @see AdvectionField_decentred_rotation
 */
template <class Mapping>
class DecentredRotationSimulation
    : public AdvectionSimulation<
              AdvectionField_decentred_rotation,
              FunctionToBeAdvected_cos_4_elipse<Mapping>>
{
public:
    /**
     * @brief Instantiate a DecentredRotationSimulation simulation.
     *
     * @param[in] mapping
     *      The mapping from the logical domain to the physical domain.
     */
    DecentredRotationSimulation(Mapping const& mapping)
        : AdvectionSimulation<
                AdvectionField_decentred_rotation,
                FunctionToBeAdvected_cos_4_elipse<Mapping>>(
                AdvectionField_decentred_rotation(),
                FunctionToBeAdvected_cos_4_elipse<Mapping>(mapping)) {};
};
