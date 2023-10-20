#pragma once

#include <geometry.hpp>

#include "irighthandside.hpp"

class KrookSourceAdaptive : public IRightHandSide
{
private:
    RhsType m_type;
    double m_extent;
    double m_stiffness;
    double m_amplitude;
    double m_density;
    double m_temperature;
    DFieldX m_mask;
    DFieldVx m_ftarget;

public:
    KrookSourceAdaptive(
            IDomainX const& gridx,
            IDomainVx const& gridvx,
            RhsType const type,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    KrookSourceAdaptive(KrookSourceAdaptive&&) = default;

    ~KrookSourceAdaptive() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

private:
    void get_amplitudes(DSpanSp amplitudes, DViewSpVx allfdistribu) const;
    void get_derivative(DSpanSpXVx df, DViewSpXVx f, DViewSpXVx f0) const;
};
