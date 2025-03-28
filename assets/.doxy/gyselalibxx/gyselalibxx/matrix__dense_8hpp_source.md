

# File matrix\_dense.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_dense.hpp**](matrix__dense_8hpp.md)

[Go to the documentation of this file](matrix__dense_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix.hpp"

class Matrix_Dense : public Matrix
{
public:
    explicit Matrix_Dense(int n);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double aij) override;

private:
    virtual int factorise_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    std::unique_ptr<int[]> ipiv;
    std::unique_ptr<double[]> a;
};
```


