

# File matrix\_corner\_block.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_corner\_block.hpp**](matrix__corner__block_8hpp.md)

[Go to the documentation of this file](matrix__corner__block_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"
#include "matrix_dense.hpp"
#include "view.hpp"

class Matrix_Corner_Block : public Matrix
{
public:
    Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual void factorise() override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const override;

protected:
    Matrix_Corner_Block(
            int n,
            int k,
            std::unique_ptr<Matrix> q,
            int lambda_size1,
            int lambda_size2);

    virtual void calculate_delta_to_factorise();

    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const;

    virtual DSpan1D solve_lambda_section_transpose(DSpan1D u, DView1D v) const;

    virtual DSpan1D solve_gamma_section(DSpan1D const u, DView1D const v) const;

    virtual DSpan1D solve_gamma_section_transpose(DSpan1D const v, DView1D const u) const;

    int const k;
    int const nb;
    std::unique_ptr<double[]> Abm_1_gamma_ptr;
    std::unique_ptr<double[]> lambda_ptr;

    std::unique_ptr<Matrix> q_block;

    Matrix_Dense delta;
    DSpan2D Abm_1_gamma;

    DSpan2D lambda;

private:
    virtual int factorise_method() override
    {
        return 0;
    }
    virtual int solve_inplace_method(double*, char, int) const override
    {
        return 0;
    }
};
```


