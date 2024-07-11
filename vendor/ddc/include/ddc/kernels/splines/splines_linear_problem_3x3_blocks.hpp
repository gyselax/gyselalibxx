// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <memory>
#include <string>

#include "splines_linear_problem.hpp"
#include "splines_linear_problem_2x2_blocks.hpp"

namespace ddc::detail {

/**
 * @brief A 3x3-blocks linear problem dedicated to the computation of a spline approximation,
 * with all blocks except center one being stored in dense format.
 *
 * A = | a | b | c |
 *     | d | Q | e |
 *     | f | g | h |
 *
 * The storage format is dense for all blocks except center one, whose storage format is determined by its type.
 *
 * The matrix itself and blocks a, Q and h are square (which fully determines the dimensions of the others).
 *
 * This class implements row & columns interchanges of the matrix and of multiple right-hand sides to restructure the
 * 3x3-blocks linear problem into a 2x2-blocks linear problem, relying then on the operations implemented in SplinesLinearProblem2x2Blocks.
 *
 * @tparam ExecSpace The Kokkos::ExecutionSpace on which operations related to the matrix are supposed to be performed.
 */
template <class ExecSpace>
class SplinesLinearProblem3x3Blocks : public SplinesLinearProblem2x2Blocks<ExecSpace>
{
public:
    using typename SplinesLinearProblem2x2Blocks<ExecSpace>::MultiRHS;
    using SplinesLinearProblem2x2Blocks<ExecSpace>::size;
    using SplinesLinearProblem2x2Blocks<ExecSpace>::solve;
    using SplinesLinearProblem2x2Blocks<ExecSpace>::m_top_left_block;

protected:
    std::size_t m_top_size;

public:
    /**
     * @brief SplinesLinearProblem3x3Blocks constructor.
     *
     * @param mat_size The size of one of the dimensions of the square matrix.
     * @param top_size The size of one of the dimensions of the top-left square block.
     * @param center_block A pointer toward the center SplinesLinearProblem. `setup_solver` must not have been called on it.
     */
    explicit SplinesLinearProblem3x3Blocks(
            std::size_t const mat_size,
            std::size_t const top_size,
            std::unique_ptr<SplinesLinearProblem<ExecSpace>> center_block)
        : SplinesLinearProblem2x2Blocks<ExecSpace>(mat_size, std::move(center_block))
        , m_top_size(top_size)
    {
    }

private:
    /// @brief Adjust indices, governs the row & columns interchanges to restructure the 3x3-blocks matrix into a 2x2-blocks matrix.
    void adjust_indices(std::size_t& i, std::size_t& j) const
    {
        std::size_t const nq = m_top_left_block->size(); // size of the center block

        if (i < m_top_size) {
            i += nq;
        } else if (i < m_top_size + nq) {
            i -= m_top_size;
        }

        if (j < m_top_size) {
            j += nq;
        } else if (j < m_top_size + nq) {
            j -= m_top_size;
        }
    }

public:
    double get_element(std::size_t i, std::size_t j) const override
    {
        adjust_indices(i, j);
        return SplinesLinearProblem2x2Blocks<ExecSpace>::get_element(i, j);
    }

    void set_element(std::size_t i, std::size_t j, double const aij) override
    {
        adjust_indices(i, j);
        return SplinesLinearProblem2x2Blocks<ExecSpace>::set_element(i, j, aij);
    }

private:
    /**
     * @brief Perform row interchanges on multiple right-hand sides to get a 2-blocks structure (matching the requirements
     * of the SplinesLinearProblem2x2Blocks solver).
     *
     * |  b_top   |    | b_center |
     * | b_center | -> |  b_top   | -- Considered as a
     * | b_bottom |    | b_bottom | -- single bottom block
     *
     * @param b The multiple right-hand sides.
     */
    void interchange_rows_from_3_to_2_blocks_rhs(MultiRHS const b) const
    {
        std::size_t const nq = m_top_left_block->size(); // size of the center block

        // prevent Kokkos::deep_copy(b_top_dst, b_top) to be a deep_copy between overlapping allocations
        assert(nq >= m_top_size);

        MultiRHS const b_top = Kokkos::
                subview(b, std::pair<std::size_t, std::size_t> {0, m_top_size}, Kokkos::ALL);
        MultiRHS const b_center = Kokkos::
                subview(b,
                        std::pair<std::size_t, std::size_t> {m_top_size, m_top_size + nq},
                        Kokkos::ALL);

        MultiRHS const b_center_dst
                = Kokkos::subview(b, std::pair<std::size_t, std::size_t> {0, nq}, Kokkos::ALL);
        MultiRHS const b_top_dst = Kokkos::
                subview(b, std::pair<std::size_t, std::size_t> {nq, nq + m_top_size}, Kokkos::ALL);

        MultiRHS const buffer = Kokkos::create_mirror(ExecSpace(), b_center);

        Kokkos::deep_copy(buffer, b_center);
        Kokkos::deep_copy(b_top_dst, b_top);
        Kokkos::deep_copy(b_center_dst, buffer);
    }

    /**
     * @brief Perform row interchanges on multiple right-hand sides to restore its 3-blocks structure.
     *
     * | b_center |    |  b_top   |
     * |  b_top   | -> | b_center |
     * | b_bottom |    | b_bottom |
     *
     * @param b The multiple right-hand sides.
     */
    void interchange_rows_from_2_to_3_blocks_rhs(MultiRHS const b) const
    {
        std::size_t const nq = m_top_left_block->size(); // size of the center block

        // prevent Kokkos::deep_copy(b_top, b_top_src) to be a deep_copy between overlapping allocations
        assert(nq >= m_top_size);

        MultiRHS const b_center_src
                = Kokkos::subview(b, std::pair<std::size_t, std::size_t> {0, nq}, Kokkos::ALL);
        MultiRHS const b_top_src = Kokkos::
                subview(b, std::pair<std::size_t, std::size_t> {nq, nq + m_top_size}, Kokkos::ALL);

        MultiRHS const b_top = Kokkos::
                subview(b, std::pair<std::size_t, std::size_t> {0, m_top_size}, Kokkos::ALL);
        MultiRHS const b_center = Kokkos::
                subview(b,
                        std::pair<std::size_t, std::size_t> {m_top_size, m_top_size + nq},
                        Kokkos::ALL);

        MultiRHS const buffer = Kokkos::create_mirror(ExecSpace(), b_center);

        Kokkos::deep_copy(buffer, b_center_src);
        Kokkos::deep_copy(b_top, b_top_src);
        Kokkos::deep_copy(b_center, buffer);
    }

public:
    /**
     * @brief Solve the multiple right-hand sides linear problem Ax=b or its transposed version A^tx=b inplace.
     *
     * Perform row interchanges on multiple right-hand sides to obtain a 2x2-blocks linear problem and call the SplinesLinearProblem2x2Blocks solver.
     *
     * @param[in, out] b A 2D Kokkos::View storing the multiple right-hand sides of the problem and receiving the corresponding solution.
     * @param transpose Choose between the direct or transposed version of the linear problem.
     */
    void solve(MultiRHS const b, bool const transpose) const override
    {
        assert(b.extent(0) == size());

        interchange_rows_from_3_to_2_blocks_rhs(b);
        SplinesLinearProblem2x2Blocks<ExecSpace>::solve(b, transpose);
        interchange_rows_from_2_to_3_blocks_rhs(b);
    }
};

} // namespace ddc::detail
