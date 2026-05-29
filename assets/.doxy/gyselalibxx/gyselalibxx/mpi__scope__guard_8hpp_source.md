

# File mpi\_scope\_guard.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpi\_scope\_guard.hpp**](mpi__scope__guard_8hpp.md)

[Go to the documentation of this file](mpi__scope__guard_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

class MpiScopeGuard
{
public:
    MpiScopeGuard(int& argc, char**& argv) noexcept;

    MpiScopeGuard(MpiScopeGuard const& rhs) = delete;

    MpiScopeGuard(MpiScopeGuard&& rhs) noexcept = delete;

    ~MpiScopeGuard() noexcept;

    MpiScopeGuard& operator=(MpiScopeGuard const& rhs) = delete;

    MpiScopeGuard& operator=(MpiScopeGuard&& rhs) noexcept = delete;
};
```


