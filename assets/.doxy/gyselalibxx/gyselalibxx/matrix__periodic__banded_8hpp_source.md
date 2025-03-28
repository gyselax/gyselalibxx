

# File matrix\_periodic\_banded.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_periodic\_banded.hpp**](matrix__periodic__banded_8hpp.md)

[Go to the documentation of this file](matrix__periodic__banded_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix_corner_block.hpp"
#include "view.hpp"

class Matrix;

class Matrix_Periodic_Banded : public Matrix_Corner_Block
{
public:
    Matrix_Periodic_Banded(int n, int kl, int ku, std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual void calculate_delta_to_factorise() override;
    virtual DSpan1D solve_lambda_section(DSpan1D v, DView1D u) const override;
    virtual DSpan1D solve_lambda_section_transpose(DSpan1D u, DView1D v) const override;
    int const kl;
    int const ku;
};
```


