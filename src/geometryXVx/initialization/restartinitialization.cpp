// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "restartinitialization.hpp"

RestartInitialization::RestartInitialization(int iter_start, double& time_start)
    : m_iter_start(iter_start)
    , m_time_start(time_start)
{
}

DSpanSpXVx RestartInitialization::operator()(DSpanSpXVx const allfdistribu) const
{
    ddc::PdiEvent("restart").with("time_saved", m_time_start).with("fdistribu", allfdistribu);
    return allfdistribu;
}