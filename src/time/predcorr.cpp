#include <cmath>
#include <iostream>

#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/TaggedVector>
#include <ddc/pdi.hpp>

#include <ipoissonsolver.hpp>
#include <ivlasovsolver.h>

#include "predcorr.h"

PredCorr::PredCorr(
        const IVlasovSolver& vlasov_solver,
        const IPoissonSolver& poisson_solver,
        double dt)
    : m_vlasov_solver(vlasov_solver)
    , m_poisson_solver(poisson_solver)
    , m_dt(dt)
{
}

void PredCorr::operator()(DSpanSpXVx fdistribu, int steps) const
{
    // efield only depends on DX
    DBlockX electric_potential(fdistribu.domain<MeshX>());

    // a 2D block of the same size as fdistribu
    DBlockSpXVx fdistribu_half_t(fdistribu.domain());

    m_poisson_solver(electric_potential, fdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;

        // computation of Electric_Potential(tn)
        m_poisson_solver(electric_potential, fdistribu);

        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("time_saved", iter_time)
                .and_with("fdistribu", fdistribu)
                .and_with("electric_potential", electric_potential);

        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu);

        // predictor
        m_vlasov_solver(fdistribu_half_t, electric_potential, m_dt / 2);

        // computation of Electric_Potential(tn+1/2)
        m_poisson_solver(electric_potential, fdistribu_half_t);

        // correction on a dt
        m_vlasov_solver(fdistribu, electric_potential, m_dt);
    }

    double const final_time = iter * m_dt;
    m_poisson_solver(electric_potential, fdistribu);
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("time_saved", final_time)
            .and_with("fdistribu", fdistribu)
            .and_with("electric_potential", electric_potential);
}
