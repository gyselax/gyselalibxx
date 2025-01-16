// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix_corner_block.hpp"
#include "view.hpp"

class Matrix;

/**
 * @brief A class representing a periodic banded matrix.
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
 */
class Matrix_Periodic_Banded : public Matrix_Corner_Block
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n matrix.
     * @param kl The number of lower diagonals.
     * @param ku The number of upper diagonals.
     * @param[in] q The sub-matrix Q describing the banded sub-matrix.
     */
    Matrix_Periodic_Banded(int n, int kl, int ku, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual void calculate_delta_to_factorize() override;
    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const override;
    virtual DSpan1D solve_lambda_section_transpose(DSpan1D u, DView1D v) const override;
    /// Number of subdiagonals.
    int const kl;
    /// Number of superdiagonals.
    int const ku;
};
