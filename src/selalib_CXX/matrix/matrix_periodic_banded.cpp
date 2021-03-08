#include "matrix_periodic_banded.h"
#include "matrix_banded.h"
#include <cassert>
#include <string.h> // for memset

Matrix_Periodic_Banded::Matrix_Periodic_Banded(int n, int kl, int ku, std::unique_ptr<Matrix> q)
    : Matrix_Corner_Block(n,max(kl,ku), std::move(q), true),
    kl(kl), ku(ku)
{
    allocate_lambda();
}

void Matrix_Periodic_Banded::allocate_lambda()
{
    lambda_ptr = new double[k*(k+1)];
    memset(lambda_ptr, 0, sizeof(double)*k*(k+1));
    lambda = new mdspan_2d(lambda_ptr, k+1, k);
}

Matrix_Periodic_Banded::~Matrix_Periodic_Banded()
{}

double Matrix_Periodic_Banded::get_element(int i, int j) const
{
    assert( i >= 0 );
    assert( i  < n );
    assert( j >= 0 );
    assert( i  < n );
    if ( i < nb && j >= nb )
    {
        int d (i-j);
        if (d >  n/2) d -= n;
        if (d < -n/2) d += n;

        if ( d < - kl || d > ku ) return 0.0;

        j -= nb;
        if ( d > 0) return (*lambda)(i,j);
        else return (*lambda)(i-nb+k+1,j);
    }
    else
    {
        return Matrix_Corner_Block::get_element(i,j);
    }
}

void Matrix_Periodic_Banded::set_element(int i, int j, double a_ij)
{
    assert( i >= 0 );
    assert( i  < n );
    assert( j >= 0 );
    assert( i  < n );
    if ( i < nb && j >= nb )
    {
        int d (i-j);
        if (d >  n/2) d -= n;
        if (d < -n/2) d += n;

        if ( d < - kl || d > ku )
        {
            assert(fabs(a_ij) < 1e-20);
            return;
        }
        j -= nb;
        if ( d > 0)
        {
            (*lambda)(i,j) = a_ij;
        }
        else
        {
            (*lambda)(i-nb+k+1,j) = a_ij;
        }
    }
    else
    {
        Matrix_Corner_Block::set_element(i, j, a_ij);
    }
}

void Matrix_Periodic_Banded::calculate_delta_to_factorize()
{
    for (int i(0); i<k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int l(0); l<=i; ++l)
        {
            double lambda_il = (*lambda)(l,i);
            for (int j(0); j < k; ++j)
            {
                double new_val(delta.get_element(j,i)-lambda_il*Abm_1_gamma(j,l));
                delta.set_element(j,i,new_val);
            }
        }
        // Lower diagonals in (*lambda)
        for (int l(i+1); l<k+1; ++l)
        {
            double lambda_il = (*lambda)(l,i);
            for (int j(0); j < k; ++j)
            {
                double new_val(delta.get_element(j,i)-lambda_il*Abm_1_gamma(j,nb-k-1+l));
                delta.set_element(j,i,new_val);
            }
        }
    }
}

void Matrix_Periodic_Banded::solve_lambda_section(mdspan_1d& u, mdspan_1d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<=i; ++j)
        {
            v(i) -= (*lambda)(j,i)*u(j);
        }
        // Lower diagonals in (*lambda)
        for (int j(i+1); j<k+1; ++j)
        {
            v(i) -= (*lambda)(j, i) * u(nb-1-k+j);
        }
    }
}

void Matrix_Periodic_Banded::solve_lambda_section(mdspan_2d& u, mdspan_2d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<=i; ++j)
        {
            for (int col(0); col < v.extent(1); ++col)
            {
                v(i,col) -= (*lambda)(j,i)*u(j,col);
            }
        }
        // Lower diagonals in (*lambda)
        for (int j(i+1); j<k+1; ++j)
        {
            for (int col(0); col < v.extent(1); ++col)
            {
                v(i,col) -= (*lambda)(j, i) * u(nb-1-k+j,col);
            }
        }
    }
}

void Matrix_Periodic_Banded::solve_lambda_section_transpose(mdspan_1d& u, mdspan_1d const& v) const
{
    for (int i(0); i< k; ++i)
    {
        // Upper diagonals in (*lambda)
        for (int j(0); j<=i; ++j)
        {
            u(j) -= (*lambda)(j,i)*v(i);
        }
        // Lower diagonals in (*lambda)
        for (int j(i+1); j<k+1; ++j)
        {
            v(nb-1-k+j) -= (*lambda)(j, i) * v(i);
        }
    }
}
