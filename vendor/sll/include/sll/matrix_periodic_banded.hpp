#ifndef MATRIX_PERIODIC_BANDED_H
#define MATRIX_PERIODIC_BANDED_H
#include <memory>

#include "sll/matrix_corner_block.hpp"
#include "sll/view.hpp"

class Matrix;

class Matrix_Periodic_Banded : public Matrix_Corner_Block
{
public:
    Matrix_Periodic_Banded(int n, int kl, int ku, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual void calculate_delta_to_factorize() override;
    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const override;
    virtual DSpan2D solve_lambda_section(DSpan2D v, DView2D u) const override;
    virtual void solve_lambda_section_transpose(DSpan1D u, DSpan1D v) const override;
    int const kl; // no. of subdiagonals
    int const ku; // no. of superdiagonals
};

#endif // MATRIX_PERIODIC_BANDED_H
