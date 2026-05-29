

# File mpi\_scope\_guard.cpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpi\_scope\_guard.cpp**](mpi__scope__guard_8cpp.md)

[Go to the documentation of this file](mpi__scope__guard_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <mpi.h>

#include "mpi_scope_guard.hpp"

MpiScopeGuard::MpiScopeGuard(int& argc, char**& argv) noexcept
{
    MPI_Init(&argc, &argv);
}

MpiScopeGuard::~MpiScopeGuard() noexcept
{
    MPI_Finalize();
}
```


