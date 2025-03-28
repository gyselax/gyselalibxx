

# File matrix\_batch\_tridiag.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_batch\_tridiag.hpp**](matrix__batch__tridiag_8hpp.md)

[Go to the documentation of this file](matrix__batch__tridiag_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
// Copyright (C) The DDC development team, see COPYRIGHT.md file
//

#pragma once

#include <Kokkos_Core.hpp>

#include "matrix_batch.hpp"

template <class ExecSpace>
class MatrixBatchTridiag : public MatrixBatch<ExecSpace>
{
public:
    using typename MatrixBatch<ExecSpace>::BatchedRHS;
    using MatrixBatch<ExecSpace>::size;
    using MatrixBatch<ExecSpace>::batch_size;

private:
    using DKokkosView2D
            = Kokkos::View<double**, Kokkos::LayoutRight, typename ExecSpace::memory_space>;
    DKokkosView2D m_subdiag;
    DKokkosView2D m_diag;
    DKokkosView2D m_uppdiag;

public:
    explicit MatrixBatchTridiag(
            const int batch_size,
            const int mat_size,
            DKokkosView2D const aa,
            DKokkosView2D const bb,
            DKokkosView2D const cc)

        : MatrixBatch<ExecSpace>(batch_size, mat_size)
        , m_subdiag(aa)
        , m_diag(bb)
        , m_uppdiag(cc)
    {
    }

    bool check_stability() const
    {
        assert(batch_size() == m_subdiag.extent(0));
        assert(size() == m_subdiag.extent(1));
        assert(batch_size() == m_diag.extent(0));
        assert(size() == m_diag.extent(1));
        assert(batch_size() == m_uppdiag.extent(0));
        assert(size() == m_uppdiag.extent(1));

        int const tmp_batch_size = batch_size();
        int const tmp_mat_size = size();
        bool is_diagdom = false;
        bool is_symmetric = false;

        DKokkosView2D subdiag_proxy = m_subdiag;
        DKokkosView2D diag_proxy = m_diag;
        DKokkosView2D uppdiag_proxy = m_uppdiag;

        Kokkos::MDRangePolicy<Kokkos::Rank<2>> batch_policy({0, 0}, {tmp_batch_size, tmp_mat_size});
        Kokkos::parallel_reduce(
                "DiagDominant",
                batch_policy,
                KOKKOS_LAMBDA(int batch_idx, int k, bool& check_diag_dom) {
                    check_diag_dom = check_diag_dom
                                     && (Kokkos::abs(subdiag_proxy(batch_idx, k))
                                                 + Kokkos::abs(uppdiag_proxy(batch_idx, k))
                                         <= Kokkos::abs(diag_proxy(batch_idx, k)));
                },
                Kokkos::LAnd<bool>(is_diagdom));
        Kokkos::parallel_reduce(
                "Symmetric",
                batch_policy,
                KOKKOS_LAMBDA(int batch_idx, int k, bool& check_sym) {
                    check_sym
                            = check_sym
                              && (Kokkos::abs(
                                          subdiag_proxy(batch_idx, k) - uppdiag_proxy(batch_idx, k))
                                  < 1e-16);
                },
                Kokkos::LAnd<bool>(is_symmetric));

        return (is_diagdom || is_symmetric);
    }

    void setup_solver() final
    {
        assert(check_stability());
    }

    void solve(BatchedRHS const b) const final
    {
        assert(batch_size() == b.extent(0));
        assert(size() == b.extent(1));

        int const tmp_batch_size = m_subdiag.extent(0);
        int const tmp_mat_size = m_subdiag.extent(1);

        Kokkos::View<
                double**,
                Kokkos::LayoutRight,
                Kokkos::DefaultExecutionSpace::memory_space> const
                cprim("cprim", batch_size(), size());
        Kokkos::View<
                double**,
                Kokkos::LayoutRight,
                Kokkos::DefaultExecutionSpace::memory_space> const
                dprim("dprim", batch_size(), size());

        DKokkosView2D subdiag_proxy = m_subdiag;
        DKokkosView2D diag_proxy = m_diag;
        DKokkosView2D uppdiag_proxy = m_uppdiag;

        Kokkos::parallel_for(
                "Tridiagonal solver",
                Kokkos::RangePolicy<ExecSpace>(0, tmp_batch_size),
                KOKKOS_LAMBDA(const int batch_idx) {
                    //ForwardStep
                    cprim(batch_idx, 0) = uppdiag_proxy(batch_idx, 0) / diag_proxy(batch_idx, 0);
                    dprim(batch_idx, 0) = b(batch_idx, 0) / diag_proxy(batch_idx, 0);
                    for (int i = 1; i < tmp_mat_size; i++) {
                        cprim(batch_idx, i)
                                = uppdiag_proxy(batch_idx, i)
                                  / (diag_proxy(batch_idx, i)
                                     - subdiag_proxy(batch_idx, i) * cprim(batch_idx, i - 1));
                        dprim(batch_idx, i)
                                = (b(batch_idx, i)
                                   - subdiag_proxy(batch_idx, i) * dprim(batch_idx, i - 1))
                                  / (diag_proxy(batch_idx, i)
                                     - subdiag_proxy(batch_idx, i) * cprim(batch_idx, i - 1));
                    }
                    //BackwardStep
                    b(batch_idx, tmp_mat_size - 1) = dprim(batch_idx, tmp_mat_size - 1);
                    for (int i = tmp_mat_size - 2; i >= 0; i--) {
                        b(batch_idx, i)
                                = dprim(batch_idx, i) - cprim(batch_idx, i) * b(batch_idx, i + 1);
                    }
                });
    }
};
```


