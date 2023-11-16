// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "restartinitialization.hpp"

RestartInitialization::RestartInitialization(int iter_start, double& time_start)
    : m_iter_start(iter_start)
    , m_time_start(time_start)
{
}

device_t<DSpanSpXVx> RestartInitialization::operator()(
        device_t<DSpanSpXVx> const allfdistribu_device) const
{
    auto allfdistribu = ddc::create_mirror_view_and_copy(allfdistribu_device.span_view());
    ddc::PdiEvent("restart").with("time_saved", m_time_start).with("fdistribu", allfdistribu);
    ddc::deepcopy(allfdistribu_device, allfdistribu);
    return allfdistribu_device;
}