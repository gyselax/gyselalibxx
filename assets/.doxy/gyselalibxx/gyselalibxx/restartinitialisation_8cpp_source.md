

# File restartinitialisation.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**restartinitialisation.cpp**](restartinitialisation_8cpp.md)

[Go to the documentation of this file](restartinitialisation_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "restartinitialisation.hpp"

RestartInitialisation::RestartInitialisation(int iter_start, double& time_start)
    : m_iter_start(iter_start)
    , m_time_start(time_start)
{
}

DFieldSpXVx RestartInitialisation::operator()(DFieldSpXVx const allfdistribu) const
{
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));
    ddc::PdiEvent("restart").with("time_saved", m_time_start).with("fdistribu", allfdistribu_host);
    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    return allfdistribu;
}
```


