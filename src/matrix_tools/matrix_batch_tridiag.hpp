// SPDX-License-Identifier: MIT
// Copyright (C) The DDC development team, see COPYRIGHT.md file
//

#pragma once

#include <Kokkos_Core.hpp>

#include "matrix_batch.hpp"

/**
 * @brief A structure for solving a set of independent tridiagonal systems using a direct method.
 * The parallelism operates on the whole collection by dispatching to threads.
 * Each problem is treated sequentially, by the tridiagonal matrix algorithm (TDMA).
 * This solver is stable for tridiagonal matrices which satisfy one of the following conditions:
 * - Diagonally Dominant.
 * - Symmetric positive-definite.
 * Diagonally Dominant property is fully checked.
 * Only symmetry property is checked, positivity-definiteness is not.
 * @tparam ExecSpace The execution space related to Kokkos.
 */
template <class ExecSpace>
class MatrixBatchTridiag : public MatrixBatch<ExecSpace>
{
public:
    using typename MatrixBatch<ExecSpace>::BatchedRHS;
    using MatrixBatch<ExecSpace>::size;
    using MatrixBatch<ExecSpace>::batch_size;

private:
    /**
     * @brief Alias for 2D double Kokkos views, LayoutRight is specified.
    */
    using DKokkosView2D
            = Kokkos::View<double**, Kokkos::LayoutRight, typename ExecSpace::memory_space>;
    DKokkosView2D m_subdiag;
    DKokkosView2D m_diag;
    DKokkosView2D m_uppdiag;

public:
    /**
     * @brief Creates an instance of the MatrixBatchTridiag class.
     * First dimension is the batch, second one refers to matrix entries indexed by line.
     * The entries aa,bb,cc are 2D Kokkos views and have the same dimensions.
     * LayoutRight: means that the "last" dimension is the contiguous one.
     * aa(batch_idx,0) and cc(batch_idx,mat_size-1) are not used for any values of batch_idx.
     *
     * @param[in] batch_size The size of the set of linear problems. 
     * @param[in] mat_size The common size of each individual matrix .
     * @param[in] aa 2d Kokkos View which stores subdiagonal components for all matrices.
     * @param[in] bb 2d Kokkos View which stores diagonal components for all matrices.
     * @param[in] cc 2d Kokkos View which stores upper diagonal components for all matrices.
     */
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

    /**
     * @brief Check if the matrices are in the stability area of the solver.
     *
     * It checks if each matrix of the batch has one of the following structures:
     * - Diagonally Dominant.
     * - Symmetric.
     *  If assertion fails, there is at least one of the matrices which does not verify these conditions. 
     *  Positivity-definiteness is not checked, it is needed to fully achieve TDMA algorithm requirement: symmetric positive definite.
     *
     * @return A boolean which indicates if the stability condition may be verified.
     */
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

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     *
     * It calls check_stability function to verify if the matrices data is in range of validity of the solver.
     *
     * The stopping criterion is a reduction factor ||Axe-b||/||b||<tol with max_iter maximum iterations.
     */
    void setup_solver() final
    {
        assert(check_stability());
    }

    /**
     * @brief Solve the batched linear problem Axe=b.
     *
     * @param[in, out] b A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions.
     */
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
