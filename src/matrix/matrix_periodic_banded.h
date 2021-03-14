#ifndef MATRIX_PERIODIC_BANDED_H
#define MATRIX_PERIODIC_BANDED_H
#include "matrix.h"
#include "matrix_corner_block.h"
#include "matrix_dense.h"

class Matrix_Periodic_Banded : public Matrix_Corner_Block
{
public:
    Matrix_Periodic_Banded(int n, int kl, int ku, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual void calculate_delta_to_factorize() override;
    virtual void solve_lambda_section(mdspan_1d& v, mdspan_1d const& u) const override;
    virtual void solve_lambda_section(mdspan_2d& v, mdspan_2d const& u) const override;
    virtual void solve_lambda_section_transpose(mdspan_1d& v, mdspan_1d const& u) const override;
    const int kl; // no. of subdiagonals
    const int ku; // no. of superdiagonals
};

#endif // MATRIX_PERIODIC_BANDED_H
