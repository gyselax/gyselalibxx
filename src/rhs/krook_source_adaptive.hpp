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

public:
    KrookSourceAdaptive(
            IDomainX const& gridx,
            IDomainVx const& gridvx,
            RhsType const type,
            RhsSolver const solver_name,
            double extent,
            double stiffness,
            double amplitude,
            double density,
            double temperature);

    KrookSourceAdaptive(KrookSourceAdaptive&&) = default;

    ~KrookSourceAdaptive() override = default;

    void rhs(DSpanVx rhs, DViewSpXVx allfdistribu, double time, IndexSpX const& ispx) const;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

private:
    double get_amplitudes(DViewSpXVx allfdistribu, IndexSpX const& ispx) const;
};
