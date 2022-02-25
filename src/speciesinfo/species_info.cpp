// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include "ddc/discretization"

#include "species_info.hpp"

using std::sqrt, std::exp;

void maxwellian_initialization(
        double const density,
        double const temperature,
        double const mean_velocity,
        DSpanVx const fMaxwellian)
{
    double const inv_sqrt_2piT = 1. / sqrt(2. * M_PI * temperature);
    IDomainVx const grid_vx = fMaxwellian.domain<IDimVx>();
    for (IndexVx const iv : grid_vx) {
        CoordVx const v = to_real(iv);
        fMaxwellian(iv) = density * inv_sqrt_2piT
                          * exp(-(v - mean_velocity) * (v - mean_velocity) / (2. * temperature));
    }
}

SpeciesInformation::SpeciesInformation(
        FieldSp<int> charge,
        FieldSp<double> mass,
        FieldSp<double> n_eq,
        FieldSp<double> T_eq,
        FieldSp<double> u_eq,
        FieldSp<double> perturb_amplitude,
        FieldSp<int> perturb_mode,
        IDomainSpXVx const& domSpXVx)
    : m_charge(std::move(charge))
    , m_mass(std::move(mass))
    , m_density_eq(std::move(n_eq))
    , m_temperature_eq(std::move(T_eq))
    , m_mean_velocity_eq(std::move(u_eq))
    , m_perturb_amplitude(std::move(perturb_amplitude))
    , m_perturb_mode(std::move(perturb_mode))
    , m_maxw_values(select<IDimSp, IDimVx>(domSpXVx))
{
    for (IndexSp const isp : get_domain<IDimSp>(n_eq)) {
        // Initialization of the Maxwellian --> fill m_maxw_values
        maxwellian_initialization(
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp),
                m_maxw_values[isp]);
    }
}
