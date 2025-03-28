

# File impitranspose.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**impitranspose.hpp**](impitranspose_8hpp.md)

[Go to the documentation of this file](impitranspose_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <mpi.h>

#include "ddc_aliases.hpp"

template <class Layout1, class Layout2>
class IMPITranspose
{
    static_assert(ddc::type_seq_same_v<
                  ddc::to_type_seq_t<typename Layout1::discrete_domain_type>,
                  ddc::to_type_seq_t<typename Layout2::discrete_domain_type>>);

public:
    using idx_range_type1 = typename Layout1::discrete_domain_type;
    using idx_range_type2 = typename Layout2::discrete_domain_type;

protected:
    MPI_Comm m_comm;

public:
    explicit IMPITranspose(MPI_Comm comm) : m_comm(comm) {}
};
```


