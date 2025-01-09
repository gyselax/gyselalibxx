// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"

/**
 * @brief A matrix class representing a banded matrix.
 */
class Matrix_Banded : public Matrix
{
public:
    /**
     * @brief A constructor for the Matrix_Banded class.
     * @param n The size of the n x n Matrix.
     * @param kl The number of lower diagonals.
     * @param ku The number of upper diagonals.
     */
    Matrix_Banded(int n, int kl, int ku);

    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual int factorize_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    /// Number of subdiagonals.
    int const kl;
    /// Number of superdiagonals.
    int const ku;
    /// Number of columns in q.
    int const c;
    /// Pivot indices
    std::unique_ptr<int[]> ipiv;
    /// banded matrix representation
    std::unique_ptr<double[]> q;
};
