#pragma once

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <rk2_solver.hpp>

class KrookSourceAdaptive : public IRightHandSide
{
private:
    RhsType m_type;
    RhsSolver m_solver_name;
    double m_extent;
    double m_stiffness;
    double m_amplitude;
    double m_density;
    double m_temperature;
    DFieldX m_mask;
    DFieldVx m_ftarget;
    std::unique_ptr<ITimeSolver> m_solver;

    double const get_amplitudes(DViewSpXVx allfdistribu, IndexSpX const ispx) const;

public:
    KrookSourceAdaptive(
            IDomainX const& gridx,
            IDomainVx const& gridvx,
            RhsType type,
            RhsSolver const solver_name,
            double const extent,
            double const stiffness,
            double const amplitude,
            double const density,
            double const temperature);

    KrookSourceAdaptive(KrookSourceAdaptive&&) = default;

    ~KrookSourceAdaptive() override = default;

    void rhs(DSpanVx rhs, DViewSpXVx allfdistribu, double const time, IndexSpX const ispx) const;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;
};
