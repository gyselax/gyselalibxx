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

void PredCorr::operator()(DSpanSpXVx fdistribu, int steps) const
{
    // efield only depends on DX
    DBlockX efield(fdistribu.domain<MeshX>());

    // a 2D block of the same size as fdistribu
    DBlockSpXVx fdistribu_half_t(fdistribu.domain());

    m_efield_solver(efield, fdistribu);

    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;

        // computation of Efield(tn)
        m_efield_solver(efield, fdistribu);

        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("iter_time", iter_time)
                .and_with("fdistribu", fdistribu)
                .and_with("efield", efield);

        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu);

        // predictor
        m_vlasov_solver(fdistribu_half_t, efield, m_dt / 2);

        // computation of Efield(tn+1/2)
        m_efield_solver(efield, fdistribu_half_t);

        // correction on a dt
        m_vlasov_solver(fdistribu, efield, m_dt);
    }

    double const final_time = iter * m_dt;
    m_efield_solver(efield, fdistribu);
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("iter_time", final_time)
            .and_with("fdistribu", fdistribu)
            .and_with("efield", efield);
}
