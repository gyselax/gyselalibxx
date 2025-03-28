

# File noperturbinitialisation.hpp

[**File List**](files.md) **>** [**geometryVparMu**](dir_9a2f28dc8f538ee0f4428810facf29b8.md) **>** [**initialisation**](dir_99d29839093a8e7b0be0d596be7efa54.md) **>** [**noperturbinitialisation.hpp**](noperturbinitialisation_8hpp.md)

[Go to the documentation of this file](noperturbinitialisation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iinitialisation.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

class NoPerturbInitialisation : public IInitialisation
{
    DConstFieldSpVparMu m_fequilibrium;

public:
    NoPerturbInitialisation(DConstFieldSpVparMu fequilibrium);

    ~NoPerturbInitialisation() override = default;

    DFieldSpVparMu operator()(DFieldSpVparMu allfdistribu) const override;
};
```


