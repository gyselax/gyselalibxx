#include "matrix_corner_block.h"
#include <cassert>
#include <string.h>

Matrix_Corner_Block::Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q)
    : Matrix(n), k(k), nb(n-k), Abm_1_gamma_ptr(new double[k*nb]),
    q_block(std::move(q)), delta(k), Abm_1_gamma(Abm_1_gamma_ptr, k, nb)
{
    assert( n  > 0 );
    assert( k >= 0 );
    assert( k <= n );
    allocate_lambda();
    memset(Abm_1_gamma_ptr, 0, sizeof(double)*k*nb);
}

Matrix_Corner_Block::Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q, bool is_virtual)
    : Matrix(n), k(k), nb(n-k), Abm_1_gamma_ptr(new double[k*nb]),
    q_block(std::move(q)), delta(k), Abm_1_gamma(Abm_1_gamma_ptr, k, nb)
{
    assert( n  > 0 );
    assert( k >= 0 );
    assert( k <= n );
    memset(Abm_1_gamma_ptr, 0, sizeof(double)*k*nb);
}

void Matrix_Corner_Block::allocate_lambda()
{
    lambda_ptr = new double[k*nb];
    memset(lambda_ptr, 0, sizeof(double)*k*nb);
    lambda = new mdspan_2d(lambda_ptr, nb, k);
}

Matrix_Corner_Block::~Matrix_Corner_Block()
{
    delete[] Abm_1_gamma_ptr;
    delete[] lambda_ptr;
    delete lambda;
}

double Matrix_Corner_Block::get_element(int i, int j) const
{
    assert( i >= 0 );
    assert( i  < n );
    assert( j >= 0 );
    assert( i  < n );
    if ( i < nb && j < nb)
    {
        return q_block->get_element(i,j);
    }
    else if ( i >= nb && j >= nb )
    {
        return delta.get_element(i-nb, j-nb);
    }
    else if ( i >= nb )
    {
        return Abm_1_gamma(i-nb,j);
    }
    else
    {
        return (*lambda)(i, j-nb);
    }
}

void Matrix_Corner_Block::set_element(int i, int j, double a_ij)
{
    assert( i >= 0 );
    assert( i  < n );
    assert( j >= 0 );
    assert( i  < n );
    if ( i < nb && j < nb)
    {
        q_block->set_element(i,j, a_ij);
    }
    else if ( i >= nb && j >= nb )
    {
        delta.set_element(i-nb, j-nb, a_ij);
    }
    else if ( i >= nb )
    {
        Abm_1_gamma(i-nb,j) = a_ij;
    }
    else
    {
        (*lambda)(i,j-nb) = a_ij;
    }
}

void Matrix_Corner_Block::calculate_delta_to_factorize()
{
    for (int i(0); i<k; ++i)
    {
        // Upper diagonals in lambda
        for (int l(0); l<nb; ++l)
        {
            double lambda_il = (*lambda)(l,i);
            for (int j(0); j < k; ++j)
            {
                double new_val(delta.get_element(j,i)-lambda_il*Abm_1_gamma(j,l));
                delta.set_element(j,i,new_val);
            }
        }
    }
}

void Matrix_Corner_Block::factorize()
{
    q_block->factorize();
    q_block->solve_inplace_matrix(Abm_1_gamma);

    calculate_delta_to_factorize();

    delta.factorize();
}

void Matrix_Corner_Block::solve_lambda_section(mdspan_1d& u, mdspan_1d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<nb; ++j)
        {
            v(i) -= (*lambda)(j,i)*u(j);
        }
    }
}
void Matrix_Corner_Block::solve_lambda_section(mdspan_2d& u, mdspan_2d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<nb; ++j)
        {
            for (int col(0); col < v.extent(1); ++col)
            {
                v(i,col) -= (*lambda)(j,i)*u(j,col);
            }
        }
    }
}

void Matrix_Corner_Block::solve_lambda_section_transpose(mdspan_1d& u, mdspan_1d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<nb; ++j)
        {
            u(j) -= (*lambda)(j,i)*v(i);
        }
    }
}

void Matrix_Corner_Block::solve_inplace(mdspan_1d& bx) const
{
    assert( bx.extent(0) == n );
    mdspan_1d u(bx.data(), nb);
    mdspan_1d v(bx.data()+nb, k);

    q_block->solve_inplace(u);

    solve_lambda_section(u, v);

    delta.solve_inplace(v);

    for (int i(0); i < nb; ++i)
    {
        double val(0);
        for (int j(0); j < k; ++j)
        {
            val += Abm_1_gamma(j,i)*v(j);
        }
        u(i)-=val;
    }
}

void Matrix_Corner_Block::solve_transpose_inplace(mdspan_1d& bx) const
{
    assert( bx.extent(0) == n );
    mdspan_1d u(bx.data(), nb);
    mdspan_1d v(bx.data()+nb, k);

    delta.solve_inplace(v);

    solve_lambda_section_transpose(u,v);

    q_block->solve_inplace(u);

    for (int j(0); j < k; ++j)
    {
        double val(0);
        for (int i(0); i < nb; ++i)
        {
            val += Abm_1_gamma(j,i)*v(i);
        }
        v(j)-=val;
    }
}

void Matrix_Corner_Block::solve_inplace_matrix(mdspan_2d& bx) const
{
    assert( bx.extent(0) == n );
    mdspan_2d u(bx.data(), nb, bx.extent(1));
    mdspan_2d v(bx.data()+nb*bx.extent(1), k, bx.extent(1));

    q_block->solve_inplace_matrix(u);

    solve_lambda_section(u, v);

    delta.solve_inplace_matrix(v);

    for (int col(0); col < bx.extent(1); ++col)
    {
        for (int i(0); i < nb; ++i)
        {
            double val(0);
            for (int j(0); j < k; ++j)
            {
                val += Abm_1_gamma(j,i)*v(j,col);
            }
            u(i,col)-=val;
        }
    }
}
