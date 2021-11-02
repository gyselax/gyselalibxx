#ifndef MATRIX_CORNER_BLOCK_H
#define MATRIX_CORNER_BLOCK_H
#include <memory>

#include "sll/matrix.hpp"
#include "sll/matrix_dense.hpp"
#include "sll/view.hpp"

class Matrix_Corner_Block : public Matrix
{
public:
    Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual void factorize() override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_inplace_matrix(DSpan2D bx) const override;

protected:
    Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q, int lambda_size);
    virtual void calculate_delta_to_factorize();
    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const;
    virtual DSpan2D solve_lambda_section(DSpan2D v, DView2D u) const;
    virtual void solve_lambda_section_transpose(DSpan1D u, DSpan1D v) const;
    int const k; // small block size
    int const nb; // main block matrix size
    std::unique_ptr<double[]> Abm_1_gamma_ptr;
    std::unique_ptr<double[]> lambda_ptr;
    //-------------------------------------
    //
    //    q = | q_block | gamma |
    //        |  lambda | delta |
    //
    //-------------------------------------
    std::unique_ptr<Matrix> q_block;
    Matrix_Dense delta;
    DSpan2D Abm_1_gamma;
    DSpan2D lambda;

private:
    virtual int factorize_method() override
    {
        return 0;
    }
    virtual int solve_inplace_method(double* b, char transpose, int nrows, int ncols) const override
    {
        return 0;
    }
};

#endif // MATRIX_CORNER_BLOCK_H
