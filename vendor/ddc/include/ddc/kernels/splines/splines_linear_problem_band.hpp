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
 * @brief A band linear problem dedicated to the computation of a spline approximation.
 *
 * The storage format is dense row-major. Lapack is used to perform every matrix and linear solver-related operations.
 *
 * Given the linear system A*x=b, we assume that A is a square (n by n)
 * with ku superdiagonals and kl subdiagonals.
 * All non-zero elements of A are stored in the rectangular matrix q, using
 * the format required by DGBTRF (LAPACK): diagonals of A are rows of q.
 * q has 2*kl rows for the subdiagonals, 1 row for the diagonal, and ku rows
 * for the superdiagonals. (The kl additional rows are needed for pivoting.)
 *
 * @tparam ExecSpace The Kokkos::ExecutionSpace on which operations related to the matrix are supposed to be performed.
 */
template <class ExecSpace>
class SplinesLinearProblemBand : public SplinesLinearProblem<ExecSpace>
{
public:
    using typename SplinesLinearProblem<ExecSpace>::MultiRHS;
    using SplinesLinearProblem<ExecSpace>::size;

protected:
    std::size_t m_kl; // no. of subdiagonals
    std::size_t m_ku; // no. of superdiagonals
    Kokkos::View<int*, Kokkos::HostSpace> m_ipiv; // pivot indices
    Kokkos::View<double**, Kokkos::HostSpace> m_q; // band matrix representation

public:
    /**
     * @brief SplinesLinearProblemBand constructor.
     *
     * @param mat_size The size of one of the dimensions of the square matrix.
     * @param kl The number of subdiagonals of the matrix.
     * @param ku The number of superdiagonals of the matrix.
     */
    explicit SplinesLinearProblemBand(
            std::size_t const mat_size,
            std::size_t const kl,
            std::size_t const ku)
        : SplinesLinearProblem<ExecSpace>(mat_size)
        , m_kl(kl)
        , m_ku(ku)
        , m_ipiv("ipiv", mat_size)
        /*
         * The matrix itself stored in band format requires a (kl + ku + 1)*mat_size 
         * allocation, but the LU-factorization requires an additional kl*mat_size block
         */
        , m_q("q", 2 * kl + ku + 1, mat_size)
    {
        assert(m_kl <= mat_size);
        assert(m_ku <= mat_size);

        Kokkos::deep_copy(m_q, 0.);
    }

private:
    std::size_t band_storage_row_index(std::size_t const i, std::size_t const j) const
    {
        return m_kl + m_ku + i - j;
    }

public:
    double get_element(std::size_t const i, std::size_t const j) const override
    {
        assert(i < size());
        assert(j < size());
        /*
         * The "row index" of the band format storage identify the (sub/super)-diagonal
         * while the column index is actually the column index of the matrix. Two layouts
         * are supported by LAPACKE. The m_kl first lines are irrelevant for the storage of 
         * the matrix itself but required for the storage of its LU factorization.
         */
        if (i >= std::
                            max(static_cast<std::ptrdiff_t>(0),
                                static_cast<std::ptrdiff_t>(j) - static_cast<std::ptrdiff_t>(m_ku))
            && i < std::min(size(), j + m_kl + 1)) {
            return m_q(band_storage_row_index(i, j), j);
        } else {
            return 0.0;
        }
    }

    void set_element(std::size_t const i, std::size_t const j, double const aij) override
    {
        assert(i < size());
        assert(j < size());
        /*
         * The "row index" of the band format storage identify the (sub/super)-diagonal
         * while the column index is actually the column index of the matrix. Two layouts
         * are supported by LAPACKE. The m_kl first lines are irrelevant for the storage of 
         * the matrix itself but required for the storage of its LU factorization.
         */
        if (i >= std::
                            max(static_cast<std::ptrdiff_t>(0),
                                static_cast<std::ptrdiff_t>(j) - static_cast<std::ptrdiff_t>(m_ku))
            && i < std::min(size(), j + m_kl + 1)) {
            m_q(band_storage_row_index(i, j), j) = aij;
        } else {
            assert(std::fabs(aij) < 1e-20);
        }
    }

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     *
     * LU-factorize the matrix A and store the pivots using the LAPACK dgbtrf() implementation.
     */
    void setup_solver() override
    {
        int const info = LAPACKE_dgbtrf(
                LAPACK_ROW_MAJOR,
                size(),
                size(),
                m_kl,
                m_ku,
                m_q.data(),
                m_q.stride(
                        0), // m_q.stride(0) if LAPACK_ROW_MAJOR, m_q.stride(1) if LAPACK_COL_MAJOR
                m_ipiv.data());
        if (info != 0) {
            throw std::runtime_error(
                    "LAPACKE_dgbtrf failed with error code " + std::to_string(info));
        }
    }

    /**
     * @brief Solve the multiple right-hand sides linear problem Ax=b or its transposed version A^tx=b inplace.
     *
     * The solver method is band gaussian elimination with partial pivoting using the LU-factorized matrix A. The implementation is LAPACK method dgbtrs.
     *
     * @param[in, out] b A 2D Kokkos::View storing the multiple right-hand sides of the problem and receiving the corresponding solution.
     * @param transpose Choose between the direct or transposed version of the linear problem.
     */
    void solve(MultiRHS b, bool const transpose) const override
    {
        assert(b.extent(0) == size());

        auto b_host = create_mirror_view(Kokkos::DefaultHostExecutionSpace(), b);
        Kokkos::deep_copy(b_host, b);
        int const info = LAPACKE_dgbtrs(
                LAPACK_ROW_MAJOR,
                transpose ? 'T' : 'N',
                b_host.extent(0),
                m_kl,
                m_ku,
                b_host.extent(1),
                m_q.data(),
                m_q.stride(
                        0), // m_q.stride(0) if LAPACK_ROW_MAJOR, m_q.stride(1) if LAPACK_COL_MAJOR
                m_ipiv.data(),
                b_host.data(),
                b_host.stride(0));
        if (info != 0) {
            throw std::runtime_error(
                    "LAPACKE_dgbtrs failed with error code " + std::to_string(info));
        }
        Kokkos::deep_copy(b, b_host);
    }
};

} // namespace ddc::detail
