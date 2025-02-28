// SPDX-License-Identifier: MIT
#pragma once

#include <iosfwd>
#include <memory>

#include "view.hpp"

/**
 * @brief The super class from which matrix classes should inherit.
 * This class is used to solve matrix equations.
 */
class Matrix
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] mat_size The size of the n x n Matrix.
     */
    Matrix(int mat_size) : n(mat_size) {}

    virtual ~Matrix() = default;

    /**
     * @brief A method to get an element from the Matrix.
     * @param[in] i The row index.
     * @param[in] j The column index.
     * @return The element at position (i,j).
     */
    virtual double get_element(int i, int j) const = 0;

    /**
     * @brief A method to set an element of the Matrix.
     * @param[in] i The row index.
     * @param[in] j The column index.
     * @param[in] a_ij The value that the element should be set to.
     */
    virtual void set_element(int i, int j, double a_ij) = 0;

    /**
     * @brief Factorise the matrix.
     * This method prepares the matrix for a call to the solve methods.
     * For most matrix types a call to factorise causes the LU decomposition
     * to be calculated. The elements of the matrix should not be accessed
     * once this method has been called.
     */
    virtual void factorise();

    /**
     * @brief Solve the matrix equation in place.
     *
     * Solve the following matrix equation:
     * @f$ M x = b @f$
     * The result @f$ x @f$ is saved into the memory allocated for @f$ b @f$.
     *
     * @param[inout] b The right-hand side of the equation on input.
     *                 The solution on output.
     * @return The solution @f$ x @f$.
     */
    virtual DSpan1D solve_inplace(DSpan1D b) const;

    /**
     * @brief Solve the transposed matrix equation in place.
     *
     * Solve the following matrix equation:
     * @f$ M^T x = b @f$
     * The result @f$ x @f$ is saved into the memory allocated for @f$ b @f$.
     *
     * @param[inout] b The right-hand side of the equation on input.
     *                 The solution on output.
     * @return The solution @f$ x @f$.
     */
    virtual DSpan1D solve_transpose_inplace(DSpan1D b) const;

    /**
     * @brief Solve multiple matrix equations in place.
     *
     * Solve the following matrix equation:
     * @f$ M x = b @f$
     * for multiple values of @f$ b @f$ and @f$ x @f$.
     * The first dimension is iterated over with each slice representing
     * an equation to be solved.
     * The result @f$ x @f$ is saved into the memory allocated for @f$ b @f$.
     *
     * @param[inout] bx The right-hand side of the equation on input.
     *                 The solution on output.
     * @return The solution @f$ x @f$.
     */
    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const;

    /**
     * @brief Get the size of the n x n Matrix.
     * @return The size of the n x n Matrix.
     */
    int get_size() const
    {
        return n;
    }

    /**
     * @brief Get a new Matrix representing a banded matrix.
     * This method returns the appropriate subclass to minimise memory and computations.
     *
     * @param[in] n The size of the n x n Matrix.
     * @param[in] kl The number of lower diagonals.
     * @param[in] ku The number of upper diagonals.
     * @param[in] pds True if the matrix is positive-definite and symmetric.
     *
     * @return An instance of a Matrix class that can store the banded matrix.
     */
    static std::unique_ptr<Matrix> make_new_banded(int n, int kl, int ku, bool pds);

    /**
     * @brief Get a new Matrix representing a periodic banded matrix.
     * This method returns the appropriate subclass to minimise memory and computations.
     * A periodic banded matrix is like a banded matrix but additionally contains non-
     * zero values in the corners. I.e it has the following sparsity pattern:
     *
     *     * * *                   * *
     *     * * * *                   *
     *     * * * * *
     *       * * * * *
     *
     *               ...
     *
     *                       * * * * *
     *     *                   * * * *
     *     * *                   * * *
     *
     * @param[in] n The size of the n x n Matrix.
     * @param[in] kl The number of lower diagonals.
     * @param[in] ku The number of upper diagonals.
     * @param[in] pds True if the matrix is positive-definite and symmetric.
     *
     * @return An instance of a Matrix class that can store the banded matrix.
     */
    static std::unique_ptr<Matrix> make_new_periodic_banded(int n, int kl, int ku, bool pds);

    /**
     * @brief Get a new Matrix representing a matrix which has a banded region.
     * This method returns the appropriate subclass to minimise memory and computations.
     * This matrix must be able to be described by the following block matrices:
     *
     *     A B C
     *     D E F
     *     G H I
     *
     * where E is a banded matrix, and A, E and I are square matrices.
     *
     * @param[in] n The size of the n x n Matrix.
     * @param[in] kl The number of lower diagonals.
     * @param[in] ku The number of upper diagonals.
     * @param[in] pds True if the matrix is positive-definite and symmetric.
     * @param[in] block1_size The size of the matrix A.
     * @param[in] block2_size The size of the matrix I.
     *
     * @return An instance of a Matrix class that can store the banded matrix.
     */
    static std::unique_ptr<Matrix> make_new_block_with_banded_region(
            int n,
            int kl,
            int ku,
            bool pds,
            int block1_size,
            int block2_size = 0);

protected:
    /**
     * @brief Call the factorisation method.
     *
     * @return The LAPACK error code.
     */
    virtual int factorise_method() = 0;

    /**
     * @brief Call the LAPACK solve method.
     * @param[inout] b The data describing the right-hand side.
     * @param[in] transpose A character ['N'/'T'] describing whether the matrix or the
     *  transposed matrix appears in the matrix equation.
     * @param[in] n_equations The number of equations being solved.
     * @return The LAPACK error code.
     */
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const = 0;

    /// The matrix size.
    int const n;
};

std::ostream& operator<<(std::ostream& o, Matrix const& m);
