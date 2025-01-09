// SPDX-License-Identifier: MIT
#pragma once

#include <memory>

#include "matrix.hpp"

/**
 * @brief A class representing a real symmetric positive definite matrix.
 */
class Matrix_PDS_Tridiag : public Matrix
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n Matrix.
     */
    explicit Matrix_PDS_Tridiag(int n);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual int factorize_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    /// The values on the diagonal.
    std::unique_ptr<double[]> d;
    /// The values on the lower diagonal.
    std::unique_ptr<double[]> l;
};
