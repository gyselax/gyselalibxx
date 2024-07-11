// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <memory>
#include <string>

#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

#include "splines_linear_problem.hpp"

namespace ddc::detail {

/**
 * @brief A positive-definite symmetric band linear problem dedicated to the computation of a spline approximation.
 *
 * The storage format is dense row-major. Lapack is used to perform every matrix and linear solver-related operations.
 *
 * Given the linear system A*x=b, we assume that A is a square (n by n)
 * with kd sub and superdiagonals.
 * All non-zero elements of A are stored in the rectangular matrix q, using
 * the format required by DPBTRF (LAPACK): (super-)diagonals of A are rows of q.
 * q has 1 row for the diagonal and kd rows for the superdiagonals.
 *
 * @tparam ExecSpace The Kokkos::ExecutionSpace on which operations related to the matrix are supposed to be performed.
 */
template <class ExecSpace>
class SplinesLinearProblemPDSBand : public SplinesLinearProblem<ExecSpace>
{
public:
    using typename SplinesLinearProblem<ExecSpace>::MultiRHS;
    using SplinesLinearProblem<ExecSpace>::size;

protected:
    Kokkos::View<double**, Kokkos::HostSpace> m_q; // pds band matrix representation

public:
    /**
     * @brief SplinesLinearProblemPDSBand constructor.
     *
     * @param mat_size The size of one of the dimensions of the square matrix.
     * @param kd The number of sub/superdiagonals of the matrix.
     */
    explicit SplinesLinearProblemPDSBand(std::size_t const mat_size, std::size_t const kd)
        : SplinesLinearProblem<ExecSpace>(mat_size)
        , m_q("q", kd + 1, mat_size)
    {
        assert(m_q.extent(0) <= mat_size);

        Kokkos::deep_copy(m_q, 0.);
    }

    double get_element(std::size_t i, std::size_t j) const override
    {
        assert(i < size());
        assert(j < size());

        // Indices are swapped for an element on subdiagonal
        if (i > j) {
            std::swap(i, j);
        }
        if (j - i < m_q.extent(0)) {
            return m_q(j - i, i);
        } else {
            return 0.0;
        }
    }

    void set_element(std::size_t i, std::size_t j, double const aij) override
    {
        assert(i < size());
        assert(j < size());

        // Indices are swapped for an element on subdiagonal
        if (i > j) {
            std::swap(i, j);
        }
        if (j - i < m_q.extent(0)) {
            m_q(j - i, i) = aij;
        } else {
            assert(std::fabs(aij) < 1e-20);
        }
    }

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     *
     * LU-factorize the matrix A and store the pivots using the LAPACK dpbtrf() implementation.
     */
    void setup_solver() override
    {
        int const info = LAPACKE_dpbtrf(
                LAPACK_ROW_MAJOR,
                'L',
                size(),
                m_q.extent(0) - 1,
                m_q.data(),
                m_q.stride(
                        0) // m_q.stride(0) if LAPACK_ROW_MAJOR, m_q.stride(1) if LAPACK_COL_MAJOR
        );
        if (info != 0) {
            throw std::runtime_error(
                    "LAPACKE_dpbtrf failed with error code " + std::to_string(info));
        }
    }

    /**
     * @brief Solve the multiple right-hand sides linear problem Ax=b or its transposed version A^tx=b inplace.
     *
     * The solver method is band gaussian elimination with partial pivoting using the LU-factorized matrix A. The implementation is LAPACK method dpbtrs.
     *
     * @param[in, out] b A 2D Kokkos::View storing the multiple right-hand sides of the problem and receiving the corresponding solution.
     * @param transpose Choose between the direct or transposed version of the linear problem (unused for a symmetric problem).
     */
    void solve(MultiRHS b, bool const) const override
    {
        assert(b.extent(0) == size());

        auto b_host = create_mirror_view(Kokkos::DefaultHostExecutionSpace(), b);
        Kokkos::deep_copy(b_host, b);
        int const info = LAPACKE_dpbtrs(
                LAPACK_ROW_MAJOR,
                'L',
                b_host.extent(0),
                m_q.extent(0) - 1,
                b_host.extent(1),
                m_q.data(),
                m_q.stride(
                        0), // m_q.stride(0) if LAPACK_ROW_MAJOR, m_q.stride(1) if LAPACK_COL_MAJOR
                b_host.data(),
                b_host.stride(0));
        if (info != 0) {
            throw std::runtime_error(
                    "LAPACKE_dpbtrs failed with error code " + std::to_string(info));
        }
        Kokkos::deep_copy(b, b_host);
    }
};

} // namespace ddc::detail
