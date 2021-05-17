#include <string.h> //for memcpy

#include "matrix_center_block.h"

Matrix_Center_Block::Matrix_Center_Block(
        int n,
        int top_block_size,
        int bottom_block_size,
        std::unique_ptr<Matrix> q)
    : Matrix_Corner_Block(n, top_block_size + bottom_block_size, std::move(q))
    , top_block_size(top_block_size)
    , bottom_block_size(bottom_block_size)
    , bottom_block_index(n - bottom_block_size)
    , swap_array(std::make_unique<double[]>(q_block->get_size()))
{
}

void Matrix_Center_Block::adjust_indexes(int& i, int& j) const
{
    if (i < top_block_size)
        i += q_block->get_size();
    else if (i < bottom_block_index)
        i -= top_block_size;

    if (j < top_block_size)
        j += q_block->get_size();
    else if (j < bottom_block_index)
        j -= top_block_size;
}

double Matrix_Center_Block::get_element(int i, int j) const
{
    adjust_indexes(i, j);
    return Matrix_Corner_Block::get_element(i, j);
}

void Matrix_Center_Block::set_element(int i, int j, double a_ij)
{
    adjust_indexes(i, j);
    Matrix_Corner_Block::set_element(i, j, a_ij);
}

void Matrix_Center_Block::swap_array_to_corner(DSpan1D& bx) const
{
    memcpy(swap_array.get(), bx.data() + top_block_size, q_block->get_size() * sizeof(double));
    memcpy(bx.data() + q_block->get_size(), bx.data(), top_block_size * sizeof(double));
    memcpy(bx.data(), swap_array.get(), q_block->get_size() * sizeof(double));
}

void Matrix_Center_Block::swap_array_to_corner(DSpan2D& bx) const
{
    int ncols(bx.extent(1));
    memcpy(swap_array.get(),
           bx.data() + top_block_size * ncols,
           q_block->get_size() * ncols * sizeof(double));
    memcpy(bx.data() + q_block->get_size() * ncols,
           bx.data(),
           top_block_size * ncols * sizeof(double));
    memcpy(bx.data(), swap_array.get(), q_block->get_size() * ncols * sizeof(double));
}

void Matrix_Center_Block::swap_array_to_center(DSpan1D& bx) const
{
    memcpy(swap_array.get(), bx.data(), q_block->get_size() * sizeof(double));
    memcpy(bx.data(), bx.data() + q_block->get_size(), top_block_size * sizeof(double));
    memcpy(bx.data() + top_block_size, swap_array.get(), q_block->get_size() * sizeof(double));
}

void Matrix_Center_Block::swap_array_to_center(DSpan2D& bx) const
{
    int ncols(bx.extent(1));
    memcpy(swap_array.get(), bx.data(), q_block->get_size() * ncols * sizeof(double));
    memcpy(bx.data(),
           bx.data() + q_block->get_size() * ncols,
           top_block_size * ncols * sizeof(double));
    memcpy(bx.data() + top_block_size * ncols,
           swap_array.get(),
           q_block->get_size() * ncols * sizeof(double));
}

void Matrix_Center_Block::solve_inplace(DSpan1D& bx) const
{
    swap_array_to_corner(bx);
    Matrix_Corner_Block::solve_inplace(bx);
    swap_array_to_center(bx);
}

void Matrix_Center_Block::solve_transpose_inplace(DSpan1D& bx) const
{
    swap_array_to_corner(bx);
    Matrix_Corner_Block::solve_transpose_inplace(bx);
    swap_array_to_center(bx);
}

void Matrix_Center_Block::solve_inplace_matrix(DSpan2D& bx) const
{
    swap_array_to_corner(bx);
    Matrix_Corner_Block::solve_inplace_matrix(bx);
    swap_array_to_center(bx);
}
