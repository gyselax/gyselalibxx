#include <cmath>
#include <iostream>

#include <ddc/ChunkSpan>
#include <ddc/DiscreteDomain>
#include <ddc/pdi.hpp>

#include <ipoissonsolver.hpp>
#include <ivlasovsolver.hpp>

#include "predcorr.hpp"

PredCorr::PredCorr(
        IVlasovSolver const& vlasov_solver,
        IPoissonSolver const& poisson_solver,
        double const dt)
    : m_vlasov_solver(vlasov_solver)
    , m_poisson_solver(poisson_solver)
    , m_dt(dt)
{
}

DSpanSpXVx PredCorr::operator()(DSpanSpXVx const allfdistribu, int const steps) const
{
    // efield only depends on DX
    DFieldX electric_potential(allfdistribu.domain<IDimX>());

    // a 2D chunck of the same size as fdistribu
    DFieldSpXVx allfdistribu_half_t(allfdistribu.domain());

    m_poisson_solver(electric_potential, allfdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;

        // computation of Electric_Potential(tn)
        m_poisson_solver(electric_potential, allfdistribu);

        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", allfdistribu)
                .and_with("electric_potential", electric_potential);

        // copy fdistribu
        deepcopy(allfdistribu_half_t, allfdistribu);

        // predictor
        m_vlasov_solver(allfdistribu_half_t, electric_potential, m_dt / 2);

        // computation of Electric_Potential(tn+1/2)
        m_poisson_solver(electric_potential, allfdistribu_half_t);

        // correction on a dt
        m_vlasov_solver(allfdistribu, electric_potential, m_dt);
    }

    double const final_time = iter * m_dt;
    m_poisson_solver(electric_potential, allfdistribu);
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", allfdistribu)
            .and_with("electric_potential", electric_potential);

    return allfdistribu;
}
