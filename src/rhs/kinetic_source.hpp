#pragma once

#include <cmath>

#include <geometry.hpp>

#include "irighthandside.hpp"

class KineticSource : public IRightHandSide
{
private:
    double m_extent;
    double m_stiffness;
    double m_amplitude;
    double m_density;
    double m_energy;
    double m_temperature;
    DFieldX m_mask_source;
    DFieldVx m_kinetic_source_vx;

public:
    KineticSource(
            IDomainX const& gridx,
            IDomainVx const& gridv,
            double const extent,
            double const stiffness,
            double const amplitude,
            double const density,
            double const energy,
            double const temperature);

    ~KineticSource() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
