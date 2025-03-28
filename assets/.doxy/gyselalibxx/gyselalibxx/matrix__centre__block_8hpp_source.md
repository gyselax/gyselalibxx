

# File matrix\_centre\_block.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_centre\_block.hpp**](matrix__centre__block_8hpp.md)

[Go to the documentation of this file](matrix__centre__block_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <memory>

#include "matrix_corner_block.hpp"
#include "view.hpp"

class Matrix;

class Matrix_Centre_Block : public Matrix_Corner_Block
{
public:
    Matrix_Centre_Block(
            int n,
            int top_block_size,
            int bottom_block_size,
            std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_multiple_inplace(DSpan2D bx) const override;

protected:
    void adjust_indexes(int& i, int& j) const;

    DSpan1D swap_array_to_corner(DSpan1D bx) const;

    DSpan1D swap_array_to_centre(DSpan1D bx) const;

    DSpan2D swap_array_to_corner(DSpan2D bx) const;

    DSpan2D swap_array_to_centre(DSpan2D bx) const;

    int const top_block_size;
    int const bottom_block_size;
    int const bottom_block_index;
    std::unique_ptr<double[]> swap_array;
};
```


