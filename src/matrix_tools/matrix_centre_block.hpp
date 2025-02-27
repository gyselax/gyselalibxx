// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix_corner_block.hpp"
#include "view.hpp"

class Matrix;

/**
 * @brief A Matrix representing a matrix which has a banded region.
 * This matrix must be able to be described by the following block matrices:
 *
 *     A B C
 *     D E F
 *     G H I
 *
 * where E is a banded matrix, and A, E and I are square matrices.
 *
 * This matrix is solved by rearranging it to:
 *
 *     E D F
 *     B A C
 *     H G I
 *
 * This new matrix is a corner matrix of the form:
 *
 *      |    Q   | gamma |
 *      | lambda | delta |
 *
 * with @f$ Q = E @f$ and
 *
 *     gamma = D F
 * 
 *     lambda = B
 *              H
 *
 *     delta = A C
 *             G I
 *
 * Internally the matrix is saved in the corner format.
 */
class Matrix_Centre_Block : public Matrix_Corner_Block
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n Matrix.
     * @param[in] top_block_size The size of the matrix A.
     * @param[in] bottom_block_size The size of the matrix I.
     * @param[in] q The banded section of the matrix.
     */
    Matrix_Centre_Block(
            int n,
            int top_block_size,
            int bottom_block_size,
            std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const override;

protected:
    /**
     * @brief Change the indices so they index the matrix in the corner layout.
     * @param[inout] i On input: the row index of the matrix with the banded region in the middle.
     *                 On output: the row index of the matrix with the banded region in the corner.
     * @param[inout] j On input: the column index of the matrix with the banded region in the middle.
     *                 On output: the column index of the matrix with the banded region in the corner.
     */
    void adjust_indexes(int& i, int& j) const;

    /**
     * @brief Rearrange the elements of the right-hand side of the equation so they are aligned
     * with the way in which the matrix is stored.
     *
     * I.e. for
     *
     *      bx = A
     *           B
     *           C
     *
     * with $f$ B @f$ a vector of length nb, return
     *
     *      bx = B
     *           A
     *           C
     *
     * @param[inout] bx The right-hand side of the matrix equation to be solved.
     * @return The right-hand side of the rearranged matrix equation.
     */
    DSpan1D swap_array_to_corner(DSpan1D bx) const;

    /**
     * @brief Rearrange the elements of the solution of the equation calculated by the Matrix_Corner_Block
     * superclass so they are aligned with the way in which the user accesses the data.
     *
     * I.e. for
     *
     *      bx = B
     *           A
     *           C
     *
     * with $f$ B @f$ a vector of length nb, return
     *
     *      bx = A
     *           B
     *           C
     *
     * @param[inout] bx The solution of the rearranged matrix equation.
     * @return The solution of the matrix equation.
     */
    DSpan1D swap_array_to_centre(DSpan1D bx) const;

    /**
     * @brief Rearrange the elements of the right-hand sides of the equations so they are aligned
     * with the way in which the matrix is stored.
     *
     * I.e. for
     *
     *      bx = A
     *           B
     *           C
     *
     * with $f$ B @f$ a vector of length nb, return
     *
     *      bx = B
     *           A
     *           C
     *
     * @param[inout] bx The right-hand sides of the matrix equations to be solved.
     * @return The right-hand sides of the rearranged matrix equations.
     */
    DSpan2D swap_array_to_corner(DSpan2D bx) const;

    /**
     * @brief Rearrange the elements of the solutions of the equations calculated by the Matrix_Corner_Block
     * superclass so they are aligned with the way in which the user accesses the data.
     *
     * I.e. for
     *
     *      bx = B
     *           A
     *           C
     *
     * with $f$ B @f$ a vector of length nb, return
     *
     *      bx = A
     *           B
     *           C
     *
     * @param[inout] bx The solutions of the rearranged matrix equations.
     * @return The solutions of the matrix equations.
     */
    DSpan2D swap_array_to_centre(DSpan2D bx) const;

    /// The size of the sub-matrix A.
    int const top_block_size;
    /// The size of the sub-matrix I
    int const bottom_block_size;
    /// The index of the first element of the sub-matrix I (size(A) + size(E)).
    int const bottom_block_index;
    /// A memory block of size top_block_size which is used during swap operations.
    std::unique_ptr<double[]> swap_array;
};
