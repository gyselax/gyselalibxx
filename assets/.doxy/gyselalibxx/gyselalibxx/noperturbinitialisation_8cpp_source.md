

# File noperturbinitialisation.cpp

[**File List**](files.md) **>** [**geometryVparMu**](dir_9a2f28dc8f538ee0f4428810facf29b8.md) **>** [**initialisation**](dir_99d29839093a8e7b0be0d596be7efa54.md) **>** [**noperturbinitialisation.cpp**](noperturbinitialisation_8cpp.md)

[Go to the documentation of this file](noperturbinitialisation_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "noperturbinitialisation.hpp"


NoPerturbInitialisation::NoPerturbInitialisation(DConstFieldSpVparMu fequilibrium)
    : m_fequilibrium(fequilibrium)
{
}

DFieldSpVparMu NoPerturbInitialisation::operator()(DFieldSpVparMu const allfdistribu) const
{
    IdxRangeSpVparMu const idxrange_spvparmu = get_idx_range(allfdistribu);

    DConstFieldSpVparMu fequilibrium_proxy = get_field(m_fequilibrium);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_spvparmu,
            KOKKOS_LAMBDA(IdxSpVparMu const ispvparmu) {
                double fdistribu_val = fequilibrium_proxy(ispvparmu);
                if (fdistribu_val < 1.e-60) {
                    fdistribu_val = 1.e-60;
                }
                allfdistribu(ispvparmu) = fdistribu_val;
            });
    return allfdistribu;
}
```


