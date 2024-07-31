// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ifluidtransportsolver.hpp"
#include "ireactionrate.hpp"

/**
 * @brief A class that solves a so-called "pressure-diffusive" fluid neutral model.
 * 
 * The equation of the model that describes the evolution of the density of 
 * neutrals can be written in dimensional units as
 *  
 * @f$\partial_t n_n + \partial_x (n_{n,eq} u_i - D_p T_n \partial_x n_n) = S_n, @f$
 * 
 * where @f$n_n(x,t)@f$ is the time and space dependent neutral density, @f$T_n@f$ is the 
 * temperature of neutrals (considered constant) and @f$u_i(x)@f$ is the ion fluid velocity. 
 * 
 * In the above equation the following definitions are used: 
 * 
 * @f$n_{n,eq} = \frac{n_i n_e K_r + n_n n_i K_{cx}}{n_i K_{cx} + n_e K_i}@f$
 * 
 * and
 * 
 * @f$D_p = \frac{1}{m_n (n_i K_{cx} + n_e K_i)}@f$
 * 
 * where @f$n_i@f$ (resp. @f$n_e@f$) is the ion (resp. electron) density and @f$m_n@f$ stands
 * for the mass of neutrals. The @f$K_i@f$, @f$K_r@f$ and @f$K_{cx}@f$ coefficients 
 * represent the reaction rates of ionization, recombination and charge-exchange reactions.
 * 
 * The density source term @f$S_n@f$ is 
 * 
 * @f$S_n = n_i n_e K_r - n_n n_e K_i.@f$
 * 
 * The pressure-diffusive equation is normalized to the relevant normalization quantities 
 * of the geometryXVx folder: 
 * - densities to a reference density @f$n_0@f$;
 * - temperatures to a reference temperature @f$T_0@f$;
 * - time normalized to the electron plasma frequency @f$\omega_{pe0} = \sqrt{n_0 e^2/(m_e \varepsilon_0)}@f$;
 * - space to the Debye length @f$\lambda_{D0} = \sqrt{\varepsilon_0 T_0 / (n_0 e^2)}@f$;
 * - ion mean velocity to the ion thermal velocity @f$v_{Ti0} = \sqrt{T_0/m_i}@f$;
 * - reaction rates to a reference rate @f$K_0@f$;
 * - masses to the electron mass @f$m_e@f$.
 * 
 * With these conventions the pressure-diffusive equation can be written with all quantities 
 * normalized as
 * 
 * @f$\partial_t n_n + \partial_x (\sqrt{\frac{m_e}{m_i}}n_{n,eq} u_i - \alpha_0 D_p T_n \partial_x n_n) = \alpha_0^{-1} S_n, @f$
 * 
 * Where @f$\alpha_0@f$ is a normalization coefficient equal to @f$\alpha_0 = \omega_{pe0}/(n_0 K_0)@f$.
 * All the terms that appear in this normalized equation keep the same expression as when writing
 * the dimensional form of the model, except that quantities involved are normalized.
 * 
 * The pressure-diffusive model is solved using a RK2 time integrator.
 * Spatial derivatives are computed using splines polynomials. 
 */
class DiffusiveNeutralSolver : public IFluidTransportSolver
{
private:
    IReactionRate const& m_charge_exchange;
    IReactionRate const& m_ionization;
    IReactionRate const& m_recombination;

    double const m_normalization_coeff;

    SplineXBuilder_1d const& m_spline_x_builder;
    SplineXEvaluator_1d const& m_spline_x_evaluator;

    DConstFieldVx const m_quadrature_coeffs;

    IdxSp find_ion(IdxRangeSp const dom_kinsp) const;

public:
    /**
     * @brief Creates an instance of the DiffusiveNeutralSolver class.
     * @param[in] charge_exchange An object that represents charge-exchange reaction rate.
     * @param[in] ionization An object that represents ionization reaction rate.
     * @param[in] recombination An object that represents recombination reaction rate.
     * @param[in] normalization_coeff A normalization coefficient for the diffusive neutral model.
     * @param[in] spline_x_builder A one-dimensional spline builder.
     * @param[in] spline_x_evaluator A one-dimensional spline evaluator.
     * @param[in] quadrature_coeffs A View referencing coefficients for a quadrature.
     */
    DiffusiveNeutralSolver(
            IReactionRate const& charge_exchange,
            IReactionRate const& ionization,
            IReactionRate const& recombination,
            double const normalization_coeff,
            SplineXBuilder_1d const& spline_x_builder,
            SplineXEvaluator_1d const& spline_x_evaluator,
            DConstFieldVx const& quadrature_coeffs);

    ~DiffusiveNeutralSolver() override = default;

    /**
     * @brief Updates the neutral fluid moments according to the pressure-diffusive neutral model. 
     * 
     * Within the pressure-diffusive model only the neutral density is evolved thus it is the 
     * only fluid moments we consider.
     *
     * @param[inout] neutrals The fluid moments describing the neutrals.
     * @param[in] allfdistribu A Span referencing a constant view of the distribution function.
     * @param[in] efield A Span referencing a constant view of the electric field.
     * @param[in] dt The time step.
     *
     * @return A span referencing the neutral fluid moments passed as argument.
     */
    DFieldSpMomX operator()(
            DFieldSpMomX neutrals,
            DConstFieldSpXVx allfdistribu,
            DConstFieldX efield,
            double dt) const override;

    /**
     * @brief Computes the expression of the time derivative of the neutral fluid moments. 
     * 
     * The expression of the time derivative is given by the equation of the pressure-diffusive 
     * neutral model, that is to say
     * 
     * @f$\partial_t n_n = - \partial_x (n_{n,eq} u_i - D_p T_n \partial_x n_n) + S_n@f$
     * 
     * This function is used by the time integrator (RK2 for instance).
     *
     * @param[inout] dn The time derivative of neutral fluid moments.
     * @param[in] n The fluid moments of the neutral species.
     * @param[in] density The plasma density (for ion and electrons).
     * @param[in] velocity The plasma mean velocity (for ion and electrons).
     * @param[in] temperature The plasma temperature (for ion and electrons).
     */
    void get_derivative(
            DFieldSpMomX dn,
            DConstFieldSpMomX n,
            DConstFieldSpX density,
            DConstFieldSpX velocity,
            DConstFieldSpX temperature) const;
};
