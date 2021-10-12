#include <cmath>
#include <iostream>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/TaggedVector>
#include <ddc/pdi.hpp>

#include <iefieldsolver.h>
#include <ivlasovsolver.h>

#include "predcorr.h"

PredCorr::PredCorr(
        const IVlasovSolver& vlasov_solver,
        const IEfieldSolver& efield_solver,
        double dt,
        double time_diag)
    : m_vlasov_solver(vlasov_solver)
    , m_efield_solver(efield_solver)
    , m_dt(dt)
    , m_time_diag(time_diag)
{
}

void PredCorr::operator()(DistributionFunction& fdistribu, double electron_mass, int steps) const
{
    // efield only depends on DX
    DBlockX efield(fdistribu.domainX());

    // a 2D block of the same size as fdistribu
    DBlockXVx fdistribu_half_t(fdistribu.domain());

    m_efield_solver(efield, fdistribu.m_values);

    double const sqrt_me_on_mspecies = sqrt(electron_mass / fdistribu.m_mass);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;

        // computation of Efield(tn)
        m_efield_solver(efield, fdistribu.m_values);

        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("iter_time", iter_time)
                .and_with("fdistribu", fdistribu.m_values)
                .and_with("efield", efield);

        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu.m_values);

        // predictor
        m_vlasov_solver(
                fdistribu_half_t,
                efield,
                fdistribu.m_charge,
                sqrt_me_on_mspecies,
                m_dt / 2);

        // computation of Efield(tn+1/2)
        m_efield_solver(efield, fdistribu_half_t);

        // correction on a dt
        m_vlasov_solver(fdistribu.m_values, efield, fdistribu.m_charge, sqrt_me_on_mspecies, m_dt);
    }

    double const final_time = iter * m_dt;
    m_efield_solver(efield, fdistribu.m_values);
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("iter_time", final_time)
            .and_with("fdistribu", fdistribu.m_values)
            .and_with("efield", efield);
}
