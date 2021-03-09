#ifndef MATRIX_CENTER_BLOCK_H
#define MATRIX_CENTER_BLOCK_H
#include "matrix_corner_block.h"

class Matrix_Center_Block : public Matrix_Corner_Block {
public:
    Matrix_Center_Block(int n,
        int top_block_size,
        int bottom_block_size,
        std::unique_ptr<Matrix> q);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;
    virtual void solve_inplace(mdspan_1d& bx) const override;
    virtual void solve_transpose_inplace(mdspan_1d& bx) const override;
    virtual void solve_inplace_matrix(mdspan_2d& bx) const override;

protected:
    void adjust_indexes(int& i, int& j) const;
    void swap_array_to_corner(mdspan_1d& bx) const;
    void swap_array_to_center(mdspan_1d& bx) const;
    void swap_array_to_corner(mdspan_2d& bx) const;
    void swap_array_to_center(mdspan_2d& bx) const;
    const int top_block_size;
    const int bottom_block_size;
    const int bottom_block_index;
    std::unique_ptr<double[]> swap_array;
};

#endif // MATRIX_CENTER_BLOCK_H
