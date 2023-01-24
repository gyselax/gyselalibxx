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
    DFieldX m_spatial_extent;
    DFieldVx m_velocity_shape;

public:
    KineticSource(
            IDomainX const& gridx,
            IDomainVx const& gridv,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double energy,
            double temperature);

    ~KineticSource() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
