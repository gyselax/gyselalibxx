#include "predcorr.h"

PredCorr::PredCorr(const IVlasovSolver& vlasov, const IEfieldSolver& efield, const MDomain1D& time)
    : m_vlasov(vlasov)
    , m_efield(efield)
    , m_time(time)
{
}


void PredCorr::operator()(DBlock2D& fdistribu, double mass_ratio) const
{
    // ex only depends on DX
    DBlock1D ex(MDomain1D {fdistribu.domain(0)});

    // a 2D block of the same size as fdistribu
    DBlock2D fdistribu_half_t(fdistribu.domain());

    for (size_t iter = m_time[0].begin(); iter < m_time[0].end(); ++iter) {
        // copy fdistribu
        fdistribu_half_t = fdistribu;

        // predictor
        m_vlasov(fdistribu_half_t, mass_ratio, m_time[0].mesh().step() / 2);

        // computation of Ex(tn+1/2)
        m_efield(ex, fdistribu_half_t);

        // correction on a dt
        m_vlasov(fdistribu, mass_ratio, m_time[0].mesh().step());

        // computation of Ex(tn+1)
        m_efield(ex, fdistribu);
    }
}
