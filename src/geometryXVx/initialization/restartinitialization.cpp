// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "restartinitialization.hpp"

RestartInitialization::RestartInitialization(int iter_start, double& time_start)
    : m_iter_start(iter_start)
    , m_time_start(time_start)
{
}

DFieldSpXVx RestartInitialization::operator()(DFieldSpXVx const allfdistribu) const
{
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));
    ddc::PdiEvent("restart").with("time_saved", m_time_start).with("fdistribu", allfdistribu_host);
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    return allfdistribu;
}
