

# File matrix.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix.hpp**](matrix_8hpp.md)

[Go to the documentation of this file](matrix_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <iosfwd>
#include <memory>

#include "view.hpp"

class Matrix
{
public:
    Matrix(int mat_size) : n(mat_size) {}

    virtual ~Matrix() = default;

    virtual double get_element(int i, int j) const = 0;

    virtual void set_element(int i, int j, double a_ij) = 0;

    virtual void factorise();

    virtual DSpan1D solve_inplace(DSpan1D b) const;

    virtual DSpan1D solve_transpose_inplace(DSpan1D b) const;

    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const;

    int get_size() const
    {
        return n;
    }

    static std::unique_ptr<Matrix> make_new_banded(int n, int kl, int ku, bool pds);

    static std::unique_ptr<Matrix> make_new_periodic_banded(int n, int kl, int ku, bool pds);

    static std::unique_ptr<Matrix> make_new_block_with_banded_region(
            int n,
            int kl,
            int ku,
            bool pds,
            int block1_size,
            int block2_size = 0);

protected:
    virtual int factorise_method() = 0;

    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const = 0;

    int const n;
};

std::ostream& operator<<(std::ostream& o, Matrix const& m);
```


