// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"
#include "matrix_dense.hpp"
#include "view.hpp"

/**
 * @brief A class representing a matrix with the following block pattern:
 *
 *      |    Q   | gamma |
 *      | lambda | delta |
 *
 * where Q is a banded matrix, and Q and delta are square matrices.
 *
 * The matrix equation is factorised using a LU decomposition. The equation
 * is then solved using this decomposition as described in section 2.5.2.1
 * of Emily Bourne's thesis [1].
 *
 * [1] Emily Bourne, "Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code". December 2022.
 */
class Matrix_Corner_Block : public Matrix
{
public:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n matrix.
     * @param[in] k The size of the k x k sub-matrix delta.
     * @param[in] q The sub-matrix Q.
     */
    Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual void factorize() override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const override;

protected:
    /**
     * @brief A constructor for the matrix.
     * @param[in] n The size of the n x n matrix.
     * @param[in] k The size of the k x k sub-matrix delta.
     * @param[in] q The sub-matrix Q.
     * @param[in] lambda_size1 The number of rows necessary to store the sub-matrix lambda.
     * @param[in] lambda_size2 The number of columns necessary to store the sub-matrix lambda.
     */
    Matrix_Corner_Block(
            int n,
            int k,
            std::unique_ptr<Matrix> q,
            int lambda_size1,
            int lambda_size2);

    /**
     * @brief Calculate the contents of the dense matrix @f$ delta' @f$ that will be factorised.
     * This is the element @f$ delta' @f$ of the blockwise LU decomposition of the matrix.
     *
     *     L = |   Q    |    0   |
     *         | lambda | delta' |
     *
     * It is defined as:
     *
     * @f$ \delta' = \delta - \lambda \beta @f$
     *
     * With @f$ \beta @f$ the solution of the equation:
     *
     * @f$ Q \beta = \gamma @f$
     *
     * which should already be stored in the variable Abm_1_gamma when this method is called.
     */
    virtual void calculate_delta_to_factorize();

    /**
     * @brief Calcuate the solution to the following equation:
     *
     * @f$ v - \lambda u @f$
     *
     * The result is saved into the argument @f$ v @f$.
     * This calculation is one of the steps necessary to solve the matrix equation.
     *
     * @param[inout] v A vector of length k describing the final values in the right-hand side of
     *              the matrix to be solved.
     * @param[in] u A vector of length (n-k). This is an intermediate result of the calculation
     *              which solves the matrix equation.
     *
     * @return The vector v solution of the equation.
     */
    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const;

    /**
     * @brief Calcuate the solution to the following equation:
     *
     * @f$ u - \lambda v @f$
     *
     * The result is saved into the argument @f$ u @f$.
     * This calculation is one of the steps necessary to solve the matrix equation.
     *
     * @param[inout] u A vector of length (n-k). This is an intermediate result of the calculation
     *              which solves the transposed matrix equation.
     * @param[inout] v A vector of length k. This is an intermediate result of the calculation
     *              which solves the transposed matrix equation.
     *
     * @return The vector u solution of the equation.
     */
    virtual DSpan1D solve_lambda_section_transpose(DSpan1D u, DView1D v) const;

    /**
     * @brief Calcuate the solution to the following equation:
     *
     * @f$ u - \beta e @f$
     *
     * The result is saved into the argument @f$ u @f$.
     * This calculation is one of the steps necessary to solve the matrix equation.
     *
     * @param[inout] u A vector of length (n-k). This is an intermediate result of the calculation
     *              which solves the matrix equation.
     * @param[inout] v A vector of length k. This is an intermediate result of the calculation
     *              which solves the matrix equation.
     *
     * @return The vector u solution of the equation.
     */
    virtual DSpan1D solve_gamma_section(DSpan1D const u, DView1D const v) const;

    /**
     * @brief Calcuate the solution to the following equation:
     *
     * @f$ v - \beta u @f$
     *
     * The result is saved into the argument @f$ v @f$.
     * This calculation is one of the steps necessary to solve the matrix equation.
     *
     * @param[inout] v A vector of length k describing the final values in the right-hand side of
     *              the matrix to be solved.
     * @param[in] u A vector of length (n-k). This is an intermediate result of the calculation
     *              which solves the transposed matrix equation.
     *
     * @return The vector v solution of the equation.
     */
    virtual DSpan1D solve_gamma_section_transpose(DSpan1D const v, DView1D const u) const;

    /// The size of the dense corner matrix delta.
    int const k;
    /// The size of the banded matrix.
    int const nb;
    /**
     * @brief Data storage for the upper-right sub-matrix.
     * This memory describes gamma before the factorize() method is called.
     * This memory describes beta after the factorize() method is called.
     */
    std::unique_ptr<double[]> Abm_1_gamma_ptr;
    /// @brief Data storage for the sub-matrix lambda.
    std::unique_ptr<double[]> lambda_ptr;

    /// The sub-matrix Q which is a banded matrix.
    std::unique_ptr<Matrix> q_block;

    /// The sub-matrix delta which is a dense matrix.
    Matrix_Dense delta;
    /**
     * @brief The upper-right sub-matrix.
     * The sub-matrix gamma before the factorize() method is called.
     * The sub-matrix beta after the factorize() method is called.
     */
    DSpan2D Abm_1_gamma;

    /// The sub-matrix lambda
    DSpan2D lambda;

private:
    virtual int factorize_method() override
    {
        return 0;
    }
    virtual int solve_inplace_method(double*, char, int) const override
    {
        return 0;
    }
};
