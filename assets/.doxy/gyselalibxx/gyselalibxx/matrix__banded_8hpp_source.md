

# File matrix\_banded.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_banded.hpp**](matrix__banded_8hpp.md)

[Go to the documentation of this file](matrix__banded_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"

class Matrix_Banded : public Matrix
{
public:
    Matrix_Banded(int n, int kl, int ku);

    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual int factorise_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    int const kl;
    int const ku;
    int const c;
    std::unique_ptr<int[]> ipiv;
    std::unique_ptr<double[]> q;
};
```


