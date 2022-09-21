#pragma once

#include <geometry.hpp>

#include "isource.hpp"

class SpeciesInformation;

class Kinetic_source : public ISource
{
private:
    double const m_px_source;
    double const m_dx_source;
    double const m_source_amplitude;
    double const m_density_amplitude;
    double const m_energy_amplitude;
    double const m_temperature_source;
    DFieldX m_mask_source;
    DFieldVx m_kinetic_source_vx;

public:
    Kinetic_source(
            IDomainX const& gridx,
            IDomainVx const& gridv,
            double const m_px_source,
            double const m_dx_source,
            double const m_source_amplitude,
            double const m_density_amplitude,
            double const m_energy_amplitude,
            double const m_temperature_source);

    ~Kinetic_source() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
