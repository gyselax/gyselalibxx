#include <ddc/BlockSpan>
#include <ddc/MDomain>
#include <ddc/ProductMDomain>
#include <ddc/TaggedVector>
#include <ddc/pdi.hpp>

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

    int iter = 0;
    for (; iter < steps; ++iter) {
        PDI_expose("iter", &iter, PDI_OUT);
        expose_to_pdi("fdistribu", fdistribu);
        expose_to_pdi("ex", ex);

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

    PDI_expose("iter", &iter, PDI_OUT);
    expose_to_pdi("fdistribu", fdistribu);
    expose_to_pdi("ex", ex);

    return fdistribu;
}
