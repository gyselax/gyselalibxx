

# File matrix\_batch.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_batch.hpp**](matrix__batch_8hpp.md)

[Go to the documentation of this file](matrix__batch_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <Kokkos_Core.hpp>

#include "view.hpp"


template <typename ExecSpace>
class MatrixBatch
{
public:
    using BatchedRHS = Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>;


private:
    std::size_t m_size;
    std::size_t m_batch_size;

protected:
    explicit MatrixBatch(const std::size_t batch_size, const std::size_t mat_size)
        : m_size(mat_size)
        , m_batch_size(batch_size)
    {
    }

public:
    virtual ~MatrixBatch() = default;

    virtual void setup_solver() = 0;

    virtual void solve(BatchedRHS b) const = 0;

    std::size_t size() const
    {
        return m_size;
    }

    std::size_t batch_size() const
    {
        return m_batch_size;
    }
};
```


