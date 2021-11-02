#ifndef MATRIX_CENTER_BLOCK_H
#define MATRIX_CENTER_BLOCK_H
#include <memory>

#include "sll/matrix_corner_block.hpp"
#include "sll/view.hpp"

class Matrix;

class Matrix_Center_Block : public Matrix_Corner_Block
{
public:
    Matrix_Center_Block(
            int n,
            int top_block_size,
            int bottom_block_size,
            std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual DSpan1D solve_inplace(DSpan1D bx) const override;
    virtual DSpan1D solve_transpose_inplace(DSpan1D bx) const override;
    virtual DSpan2D solve_inplace_matrix(DSpan2D bx) const override;

protected:
    void adjust_indexes(int& i, int& j) const;
    DSpan1D swap_array_to_corner(DSpan1D bx) const;
    DSpan1D swap_array_to_center(DSpan1D bx) const;
    DSpan2D swap_array_to_corner(DSpan2D bx) const;
    DSpan2D swap_array_to_center(DSpan2D bx) const;
    int const top_block_size;
    int const bottom_block_size;
    int const bottom_block_index;
    std::unique_ptr<double[]> swap_array;
};

#endif // MATRIX_CENTER_BLOCK_H
