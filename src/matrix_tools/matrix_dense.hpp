// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"

/**
 * @brief A class describing a dense matrix.
 */
class Matrix_Dense : public Matrix
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n Matrix.
     */
    explicit Matrix_Dense(int n);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double aij) override;

private:
    virtual int factorize_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    std::unique_ptr<int[]> ipiv;
    std::unique_ptr<double[]> a;
};
