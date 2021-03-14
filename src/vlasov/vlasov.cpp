#include "vlasov.h"

Vlasov::Vlasov(const Advection1D& advec_x, const Advection1D& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


void Vlasov::operator()(DBlock2D& data, double dt) const
{
    m_advec_x(data, dt / 2);
    m_advec_vx(data, dt);
    m_advec_x(data, dt / 2);
}
