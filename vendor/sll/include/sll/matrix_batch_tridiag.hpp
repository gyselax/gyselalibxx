// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <sll/matrix_batch.hpp>

#include <Kokkos_Core.hpp>

/**
 * @brief A structure for solving a set of independant tridiagonal systems using a direct method.
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
private:
    /**
     * @brief Alias for 2D double Kokkos views, LayoutRight is specified.
    */
    using DKokkosView2D
            = Kokkos::View<double**, Kokkos::LayoutRight, typename ExecSpace::memory_space>;
    using MatrixBatch<ExecSpace>::get_size;
    using MatrixBatch<ExecSpace>::get_batch_size;
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

    MatrixBatchTridiag(
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
        assert(get_batch_size() == m_subdiag.extent(0));
        assert(get_size() == m_subdiag.extent(1));
        assert(get_batch_size() == m_diag.extent(0));
        assert(get_size() == m_diag.extent(1));
        assert(get_batch_size() == m_uppdiag.extent(0));
        assert(get_size() == m_uppdiag.extent(1));

        int const tmp_batch_size = get_batch_size();
        int const tmp_mat_size = get_size();
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
     * @brief Solves the collection of problems to find x_i such that: A_i x_i=rhs_i  i=0,batch_size.
     * @param[inout] rhs_view 2d Kokkos View which stores right hand side, it is also
     *               used to store the results.
     * @return  The computation result, stored in a 2d Kokkos View.
     */
    DKokkosView2D solve_inplace(DKokkosView2D const rhs_view) const
    {
        assert(get_batch_size() == rhs_view.extent(0));
        assert(get_size() == rhs_view.extent(1));

        int const tmp_batch_size = m_subdiag.extent(0);
        int const tmp_mat_size = m_subdiag.extent(1);

        Kokkos::View<
                double**,
                Kokkos::LayoutRight,
                Kokkos::DefaultExecutionSpace::memory_space> const
                cprim("cprim", get_batch_size(), get_size());
        Kokkos::View<
                double**,
                Kokkos::LayoutRight,
                Kokkos::DefaultExecutionSpace::memory_space> const
                dprim("dprim", get_batch_size(), get_size());

        DKokkosView2D subdiag_proxy = m_subdiag;
        DKokkosView2D diag_proxy = m_diag;
        DKokkosView2D uppdiag_proxy = m_uppdiag;

        Kokkos::parallel_for(
                "Tridiagonal solver",
                Kokkos::RangePolicy<ExecSpace>(0, tmp_batch_size),
                KOKKOS_LAMBDA(const int batch_idx) {
                    //ForwardStep
                    cprim(batch_idx, 0) = uppdiag_proxy(batch_idx, 0) / diag_proxy(batch_idx, 0);
                    dprim(batch_idx, 0) = rhs_view(batch_idx, 0) / diag_proxy(batch_idx, 0);
                    for (int i = 1; i < tmp_mat_size; i++) {
                        cprim(batch_idx, i)
                                = uppdiag_proxy(batch_idx, i)
                                  / (diag_proxy(batch_idx, i)
                                     - subdiag_proxy(batch_idx, i) * cprim(batch_idx, i - 1));
                        dprim(batch_idx, i)
                                = (rhs_view(batch_idx, i)
                                   - subdiag_proxy(batch_idx, i) * dprim(batch_idx, i - 1))
                                  / (diag_proxy(batch_idx, i)
                                     - subdiag_proxy(batch_idx, i) * cprim(batch_idx, i - 1));
                    }
                    //BackwardStep
                    rhs_view(batch_idx, tmp_mat_size - 1) = dprim(batch_idx, tmp_mat_size - 1);
                    for (int i = tmp_mat_size - 2; i >= 0; i--) {
                        rhs_view(batch_idx, i) = dprim(batch_idx, i)
                                                 - cprim(batch_idx, i) * rhs_view(batch_idx, i + 1);
                    }
                });
        return rhs_view;
    }
    /**
     * @brief function used to execute specific implementation details.
     * It calls check_stability function to verify if the matrices data is in range of validity of the solver.
    */
    void factorize() final
    {
        assert(check_stability());
    }
};
