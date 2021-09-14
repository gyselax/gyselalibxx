#include <ddc/deepcopy.h>
#include <ddc/mdomain.h>
#include <ddc/product_mdomain.h>
#include <ddc/taggedvector.h>

#include <iefieldsolver.h>
#include <ivlasovsolver.h>

#include "predcorr.h"

PredCorr::PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, RLengthT dt)
    : m_vlasov(vlasov)
    , m_efield(efield)
    , m_dt(dt)
{
}


DSpanXVx PredCorr::operator()(DSpanXVx fdistribu, double mass_ratio, int steps) const
{
    // ex only depends on DX
    DBlockX ex(fdistribu.domain<MeshX>());

    // a 2D block of the same size as fdistribu
    DBlockXVx fdistribu_half_t(fdistribu.domain());

    for (int iter = 0; iter < steps; ++iter) {
        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu);

        // predictor
        m_vlasov(fdistribu_half_t, mass_ratio, m_dt / 2);

        // computation of Ex(tn+1/2)
        m_efield(ex, fdistribu_half_t);

        // correction on a dt
        m_vlasov(fdistribu, mass_ratio, m_dt);

        // computation of Ex(tn+1)
        m_efield(ex, fdistribu);
    }

    return fdistribu;
}
