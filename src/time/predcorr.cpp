#include <ddc/BlockSpan>
#include <ddc/ProductMDomain>
#include <ddc/TaggedVector>
#include <ddc/pdi.hpp>

#include <iefieldsolver.h>
#include <ivlasovsolver.h>

#include "predcorr.h"

PredCorr::PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, double dt)
    : m_vlasov(vlasov)
    , m_efield(efield)
    , m_dt(dt)
{
}

void PredCorr::operator()(DistributionFunction& fdistribu, double electron_mass, int steps) const
{
    // ex only depends on DX
    DBlockX ex(fdistribu.values.domain<MeshX>());

    // a 2D block of the same size as fdistribu
    DBlockXVx fdistribu_half_t(fdistribu.values.domain());

    double const sqrt_me_on_mspecies = sqrt(electron_mass / fdistribu.mass);
    int iter = 0;
    for (; iter < steps; ++iter) {
        double const iter_time = iter * m_dt;
        PdiEvent("iteration")
                .with("iter", iter)
                .and_with("iter_time", iter_time)
                .and_with("fdistribu", fdistribu.values)
                .and_with("ex", ex);

        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu.values);

        // predictor
        m_vlasov(fdistribu_half_t, ex, sqrt_me_on_mspecies, m_dt / 2);

        // computation of Ex(tn+1/2)
        m_efield(ex, fdistribu_half_t);

        // correction on a dt
        m_vlasov(fdistribu.values, ex, sqrt_me_on_mspecies, m_dt);

        // computation of Ex(tn+1)
        m_efield(ex, fdistribu.values);
    }

    double const final_time = iter * m_dt;
    PdiEvent("last_iteration")
            .with("iter", iter)
            .and_with("iter_time", final_time)
            .and_with("fdistribu", fdistribu.values)
            .and_with("ex", ex);
}
