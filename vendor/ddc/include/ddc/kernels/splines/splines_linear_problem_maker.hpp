// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <optional>

#include "splines_linear_problem_2x2_blocks.hpp"
#include "splines_linear_problem_3x3_blocks.hpp"
#include "splines_linear_problem_band.hpp"
#include "splines_linear_problem_dense.hpp"
#include "splines_linear_problem_pds_band.hpp"
#include "splines_linear_problem_pds_tridiag.hpp"
#include "splines_linear_problem_sparse.hpp"

namespace ddc::detail {
/**
 * @brief A static factory for SplinesLinearProblem classes.
 */
class SplinesLinearProblemMaker
{
public:
    /**
     * @brief Construct a dense matrix
     *
     * @tparam the Kokkos::ExecutionSpace on which matrix-related operation will be performed.
     * @param n The size of one of the dimensions of the square matrix.
     *
     * @return The SplinesLinearProblem instance.
     */
    template <typename ExecSpace>
    static std::unique_ptr<SplinesLinearProblem<ExecSpace>> make_new_dense(int const n)
    {
        return std::make_unique<SplinesLinearProblemDense<ExecSpace>>(n);
    }

    /**
     * @brief Construct a band matrix
     *
     * @tparam the Kokkos::ExecutionSpace on which matrix-related operation will be performed.
     * @param n The size of one of the dimensions of the square matrix.
     * @param kl The number of subdiagonals.
     * @param ku The number of superdiagonals.
     * @param pds A boolean indicating if the matrix is positive-definite symetric or not.
     *
     * @return The SplinesLinearProblem instance.
     */
    template <typename ExecSpace>
    static std::unique_ptr<SplinesLinearProblem<ExecSpace>> make_new_band(
            int const n,
            int const kl,
            int const ku,
            bool const pds)
    {
        if (kl == ku && kl == 1 && pds) {
            return std::make_unique<SplinesLinearProblemPDSTridiag<ExecSpace>>(n);
        } else if (kl == ku && pds) {
            return std::make_unique<SplinesLinearProblemPDSBand<ExecSpace>>(n, kl);
        } else if (2 * kl + ku + 1 >= n) {
            return std::make_unique<SplinesLinearProblemDense<ExecSpace>>(n);
        } else {
            return std::make_unique<SplinesLinearProblemBand<ExecSpace>>(n, kl, ku);
        }
    }

    /**
     * @brief Construct a 2x2-blocks or 3x3-blocks linear problem with band "main" block (the one called
     * Q in SplinesLinearProblem2x2Blocks and SplinesLinearProblem3x3Blocks).
     *
     * @tparam the Kokkos::ExecutionSpace on which matrix-related operation will be performed.
     * @param n The size of one of the dimensions of the whole square matrix.
     * @param kl The number of subdiagonals in the band block.
     * @param ku The number of superdiagonals in the band block.
     * @param pds A boolean indicating if the band block is positive-definite symetric or not.
     * @param bottom_right_size The size of one of the dimensions of the bottom-right square block.
     * @param top_left_size The size of one of the dimensions of the top-left square block.
     *
     * @return The SplinesLinearProblem instance.
     */
    template <typename ExecSpace>
    static std::unique_ptr<SplinesLinearProblem<ExecSpace>>
    make_new_block_matrix_with_band_main_block(
            int const n,
            int const kl,
            int const ku,
            bool const pds,
            int const bottom_right_size,
            int const top_left_size = 0)
    {
        int const main_size = n - top_left_size - bottom_right_size;
        std::unique_ptr<SplinesLinearProblem<ExecSpace>> main_block
                = make_new_band<ExecSpace>(main_size, kl, ku, pds);
        if (top_left_size == 0) {
            return std::make_unique<
                    SplinesLinearProblem2x2Blocks<ExecSpace>>(n, std::move(main_block));
        }
        return std::make_unique<
                SplinesLinearProblem3x3Blocks<ExecSpace>>(n, top_left_size, std::move(main_block));
    }

    /**
     * @brief Construct a 2x2-blocks linear problem with band "main" block (the one called
     * Q in SplinesLinearProblem2x2Blocks) and other blocks containing the "periodic parts" of
     * a periodic band matrix.
     *
     * It simply calls make_new_block_matrix_with_band_main_block with bottom_size being
     * max(kl, ku) (except if the allocation would be higher than instantiating a SplinesLinearProblemDense).
     *
     * @tparam the Kokkos::ExecutionSpace on which matrix-related operation will be performed.
     * @param n The size of one of the dimensions of the whole square matrix.
     * @param kl The number of subdiagonals in the band block.
     * @param ku The number of superdiagonals in the band block.
     * @param pds A boolean indicating if the band block is positive-definite symetric or not.
     *
     * @return The SplinesLinearProblem instance.
     */
    template <typename ExecSpace>
    static std::unique_ptr<SplinesLinearProblem<ExecSpace>> make_new_periodic_band_matrix(
            int const n,
            int const kl,
            int const ku,
            bool const pds)
    {
        assert(kl < n);
        assert(ku < n);
        int const bottom_size = std::max(kl, ku);
        int const top_size = n - bottom_size;

        if (bottom_size * (n + top_size) + (2 * kl + ku + 1) * top_size >= n * n) {
            return std::make_unique<SplinesLinearProblemDense<ExecSpace>>(n);
        }

        return make_new_block_matrix_with_band_main_block<ExecSpace>(n, kl, ku, pds, bottom_size);
    }

    /**
     * @brief Construct a sparse matrix
     *
     * @tparam the Kokkos::ExecutionSpace on which matrix-related operation will be performed.
     * @param n The size of one of the dimensions of the square matrix.
     * @param cols_per_chunk A parameter used by the slicer (internal to the solver) to define the size
     * of a chunk of right-hand sides of the linear problem to be computed in parallel (chunks are treated
     * by the linear solver one-after-the-other).
     * This value is optional. If no value is provided then the default value is chosen by the requested solver.
     *
     * @param preconditionner_max_block_size A parameter used by the slicer (internal to the solver) to
     * define the size of a block used by the Block-Jacobi preconditioner.
     * This value is optional. If no value is provided then the default value is chosen by the requested solver.
     *
     * @return The SplinesLinearProblem instance.
     */
    template <typename ExecSpace>
    static std::unique_ptr<SplinesLinearProblem<ExecSpace>> make_new_sparse(
            int const n,
            std::optional<std::size_t> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
    {
        return std::make_unique<SplinesLinearProblemSparse<
                ExecSpace>>(n, cols_per_chunk, preconditionner_max_block_size);
    }
};

} // namespace ddc::detail
