#include "predcorr.h"

PredCorr::PredCorr(const Vlasov& vlasov, const EfieldSolver& efield)
    : m_vlasov(vlasov)
    , m_efield(efield)
{
}


void PredCorr::operator()(DBlock2D& fdistribu, const MDomain1D& time) const
{
    // ex only depends on DX
    DBlock1D ex {fdistribu.domain().slice<0>()};

    // a 2D block of the same size as fdistribu
    DBlock2D fdistribu_half_t(fdistribu.domain());

    for (size_t iter = time.begin()[0]; iter < time.end()[0]; ++iter) {
        // copy fdistribu
        fdistribu_half_t = fdistribu;

        // predictor
        m_vlasov(fdistribu_half_t, time.mesher(0).step() / 2);

        // computation of Ex(tn+1/2)
        m_efield(ex, fdistribu_half_t);

        // correction on a dt
        m_vlasov(fdistribu, time.mesher(0).step());

        // computation of Ex(tn+1)
        m_efield(ex, fdistribu);
    }
}
