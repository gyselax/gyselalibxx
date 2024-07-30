#pragma once

#include <cmath>

#include <geometry.hpp>

#include "ikineticfluidcoupling.hpp"
#include "ireactionrate.hpp"

/**
 * @brief A class that describes a source of particles due to neutrals.
 * 
 * The KineticFluidCouplingSource class solves the following evolution equations: 
 * @f$df/dt = S_n,N(x) * S_v(x,v)@f$
 * where @f$S_n,N(x)@f$ is what we call the density_source_neutral
 * @f$dn_N/dt = - S_n,N(x)@f$
 * Where @f$S_n,N(x) = n_N(x) n_e(x) K_i(x) - n_i(x) n_e(x) K_r(x)@f$
 * @f$S_v(x,v)@f$ is the sum of order 0 to 2 Hermite polynomials times a Maxwellian velocity distribution function.
 * 
 *
 * The velocity_shape_source @f$S_v(x,v)@f$ defines the velocity profile of the source in the parallel velocity direction. 
 * It is the sum of a source that injects only density, a source that injects only momentum and a source that injects only energy. 
 * If the density and energy parameters are equal to one (usual case), the resulting velocity_shape is maxwellian. 
 * 
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/kinetic_source.pdf). 
 */
class KineticFluidCouplingSource : public IKineticFluidCoupling
{
private:
    double m_density_coupling_coeff;
    double m_momentum_coupling_coeff;
    double m_energy_coupling_coeff;
    IReactionRate const& m_ionization;
    IReactionRate const& m_recombination;
    double m_normalization_coeff;
    DViewVx const m_quadrature_coeffs;

public:
    /**
     * @brief Creates an instance of the KineticFluidCouplingSource class.
     * 
     * @param[in] density_coupling_coeff The coefficient of the density source.
     * @param[in] momentum_coupling_coeff The coefficient of the momentum source.
     * @param[in] energy_coupling_coeff The coefficient of the energy source.
     * @param[in] ionization The rate of the ionization reaction.
     * @param[in] recombination The rate of the recombination reaction.
     * @param[in] normalization_coeff The normalization coefficient of neutrals.
     * @param[in] quadrature_coeffs A View referencing coefficients for a quadrature.
     */
    KineticFluidCouplingSource(
            double density_coupling_coeff,
            double momentum_coupling_coeff,
            double energy_coupling_coeff,
            IReactionRate const& ionization,
            IReactionRate const& recombination,
            double normalization_coeff,
            DViewVx const& quadrature_coeffs);

    ~KineticFluidCouplingSource() override = default;

    /**
     * @brief Update the distribution function and neutral density following the KineticFluidCouplingSource operator. 
     *
     * @param[in, out] allfdistribu The distribution function.
     * @param[in, out] neutrals The neutral density.
     * @param[in] dt The time step over which the collisions occur.
     *
     */
    void operator()(DSpanSpXVx const allfdistribu, DSpanSpMX neutrals, double const dt)
            const override;

public:
    /**
    * @brief Returns the index of the ion species in the domain.
    *
    * @param[in] dom_kinsp The domain of the kinetic species.
    * 
    * @return The index of the ion species in the domain.
    */
    IndexSp find_ion(IDomainSp const dom_kinsp) const;

    /**
     * @brief Computes the source term density_source_neutral(x) of the KineticFluidCouplingSource operator.
     * 
     * @param[in, out] density_source_neutral The source term.
     * @param[in] kinsp_density The computed plasma densities of the distribution function.
     * @param[in] neutrals The neutral density.
     * @param[in] ionization The ionization rate.
     * @param[in] recombination The recombination rate.
     * 
    */
    void get_source_term(
            DSpanX density_source_neutral,
            DViewSpX kinsp_density,
            DViewSpMX neutrals,
            DViewSpX ionization,
            DViewSpX recombination) const;

    /**
     * @brief Computes dn (of neutral density) for the equation dn/dt = - density_source_neutral.
     * 
     * @param[in, out] dn The infinitesimal of the neutral density.
     * @param[in] neutrals The neutral density.
     * @param[in] density_source_neutral The density source term.
     * 
    */
    void get_derivative_neutrals(DSpanSpMX dn, DViewSpMX neutrals, DViewX density_source_neutral)
            const;

    /**
     * @brief Computes df for the equation df/dt = density_source_neutral(x) * velocity_shape_source(x,v).
     * 
     * @param[in, out] df The infinitesimal of the distribution function.
     * @param[in] allfdistribu The distribution function.
     * @param[in] velocity_shape_source The velocity shape of the source.
     *
    */
    void get_derivative_allfdistribu(
            DSpanSpXVx df,
            DViewSpXVx allfdistribu,
            DViewSpXVx velocity_shape_source) const;
};
