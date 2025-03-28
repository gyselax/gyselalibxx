

# File restartinitialisation.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**restartinitialisation.hpp**](restartinitialisation_8hpp.md)

[Go to the documentation of this file](restartinitialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "iinitialisation.hpp"
#include "species_info.hpp"

class RestartInitialisation : public IInitialisation
{
private:
    int m_iter_start; /* iteration number to perform the restart from */
    double& m_time_start; /* corresponding simulation time */

public:
    RestartInitialisation(int iter_start, double& time_start);

    ~RestartInitialisation() override = default;

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const override;
};
```


