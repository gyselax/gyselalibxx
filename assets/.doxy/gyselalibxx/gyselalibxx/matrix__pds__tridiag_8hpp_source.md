

# File matrix\_pds\_tridiag.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_pds\_tridiag.hpp**](matrix__pds__tridiag_8hpp.md)

[Go to the documentation of this file](matrix__pds__tridiag_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <memory>

#include "matrix.hpp"

class Matrix_PDS_Tridiag : public Matrix
{
public:
    explicit Matrix_PDS_Tridiag(int n);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual int factorise_method() override;
    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override;
    std::unique_ptr<double[]> d;
    std::unique_ptr<double[]> l;
};
```


