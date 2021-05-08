#include "predcorr2.h"

PredCorr2::PredCorr2(
        const IVlasovSolver2& vlasov,
        const IEfieldSolver2& efield,
        const MDomain<Dim::T>& time)
    : m_vlasov(vlasov)
    , m_efield(efield)
    , m_time(time)
{
}


DBlockViewXVx& PredCorr2::operator()(DBlockViewXVx& fdistribu, double mass_ratio) const
{
    // ex only depends on DX
    DBlockX ex(fdistribu.domain<Dim::X>());

    // a 2D block of the same size as fdistribu
    DBlockXVx fdistribu_half_t(fdistribu.domain());

    for (auto&& iter : m_time) {
        // copy fdistribu
        deepcopy(fdistribu_half_t, fdistribu);

        // predictor
        m_vlasov(fdistribu_half_t, mass_ratio, m_time.step() / 2);

        // computation of Ex(tn+1/2)
        m_efield(ex, fdistribu_half_t);

        // correction on a dt
        m_vlasov(fdistribu, mass_ratio, m_time.step());

        // computation of Ex(tn+1)
        m_efield(ex, fdistribu);
    }

    return fdistribu;
}
