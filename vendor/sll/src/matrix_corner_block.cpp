#include <cassert>
#include <utility>

#include <Kokkos_Core.hpp>
#include <string.h>

#include "sll/matrix_corner_block.hpp"

Matrix_Corner_Block::Matrix_Corner_Block(int const n, int const k, std::unique_ptr<Matrix> q)
    : Matrix(n)
    , k(k)
    , nb(n - k)
    , Abm_1_gamma_ptr(std::make_unique<double[]>(k * nb))
    , lambda_ptr(std::make_unique<double[]>(k * nb))
    , q_block(std::move(q))
    , delta(k)
    , Abm_1_gamma(Abm_1_gamma_ptr.get(), k, nb)
    , lambda(lambda_ptr.get(), nb, k)
{
    assert(n > 0);
    assert(k >= 0);
    assert(k <= n);
    assert(nb == q_block->get_size());
    memset(lambda_ptr.get(), 0, sizeof(double) * k * nb);
    memset(Abm_1_gamma_ptr.get(), 0, sizeof(double) * k * nb);
}

Matrix_Corner_Block::Matrix_Corner_Block(
        int const n,
        int const k,
        std::unique_ptr<Matrix> q,
        int const lambda_size1,
        int const lambda_size2)
    : Matrix(n)
    , k(k)
    , nb(n - k)
    , Abm_1_gamma_ptr(std::make_unique<double[]>(k * nb))
    , lambda_ptr(std::make_unique<double[]>(lambda_size1 * lambda_size2))
    , q_block(std::move(q))
    , delta(k)
    , Abm_1_gamma(Abm_1_gamma_ptr.get(), k, nb)
    , lambda(lambda_ptr.get(), lambda_size1, lambda_size2)
{
    assert(n > 0);
    assert(k >= 0);
    assert(k <= n);
    assert(nb == q_block->get_size());
    memset(lambda_ptr.get(), 0, sizeof(double) * lambda_size1 * lambda_size2);
    memset(Abm_1_gamma_ptr.get(), 0, sizeof(double) * k * nb);
}

double Matrix_Corner_Block::get_element(int const i, int const j) const
{
    assert(i >= 0);
    assert(i < n);
    assert(j >= 0);
    assert(i < n);
    if (i < nb && j < nb) {
        return q_block->get_element(i, j);
    } else if (i >= nb && j >= nb) {
        return delta.get_element(i - nb, j - nb);
    } else if (j >= nb) {
        return Abm_1_gamma(j - nb, i);
    } else {
        return lambda(j, i - nb);
    }
}

void Matrix_Corner_Block::set_element(int const i, int const j, double const a_ij)
{
    assert(i >= 0);
    assert(i < n);
    assert(j >= 0);
    assert(i < n);
    if (i < nb && j < nb) {
        q_block->set_element(i, j, a_ij);
    } else if (i >= nb && j >= nb) {
        delta.set_element(i - nb, j - nb, a_ij);
    } else if (j >= nb) {
        Abm_1_gamma(j - nb, i) = a_ij;
    } else {
        lambda(j, i - nb) = a_ij;
    }
}

void Matrix_Corner_Block::calculate_delta_to_factorize()
{
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            double val = 0.0;
            for (int l = 0; l < nb; ++l) {
                val += lambda(l, i) * Abm_1_gamma(j, l);
            }
            delta.set_element(i, j, delta.get_element(i, j) - val);
        }
    }
}

void Matrix_Corner_Block::factorize()
{
    q_block->factorize();
    q_block->solve_multiple_inplace(Abm_1_gamma);

    calculate_delta_to_factorize();

    delta.factorize();
}

DSpan1D Matrix_Corner_Block::solve_lambda_section(DSpan1D const v, DView1D const u) const
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in lambda
        for (int j = 0; j < nb; ++j) {
            v(i) -= lambda(j, i) * u(j);
        }
    }
    return v;
}

DSpan1D Matrix_Corner_Block::solve_lambda_section_transpose(DSpan1D const u, DView1D const v) const
{
    for (int i = 0; i < nb; ++i) {
        // Upper diagonals in lambda
        for (int j = 0; j < k; ++j) {
            u(i) -= lambda(i, j) * v(j);
        }
    }
    return u;
}

DSpan1D Matrix_Corner_Block::solve_gamma_section(DSpan1D const u, DView1D const v) const
{
    for (int i = 0; i < nb; ++i) {
        double val = 0.;
        for (int j = 0; j < k; ++j) {
            val += Abm_1_gamma(j, i) * v(j);
        }
        u(i) -= val;
    }
    return u;
}

DSpan1D Matrix_Corner_Block::solve_gamma_section_transpose(DSpan1D const v, DView1D const u) const
{
    for (int j = 0; j < k; ++j) {
        double val = 0.;
        for (int i = 0; i < nb; ++i) {
            val += Abm_1_gamma(j, i) * u(i);
        }
        v(j) -= val;
    }
    return v;
}

DSpan1D Matrix_Corner_Block::solve_inplace(DSpan1D const bx) const
{
    assert(int(bx.extent(0)) == n);

    DSpan1D const u(bx.data_handle(), nb);
    DSpan1D const v(bx.data_handle() + nb, k);

    q_block->solve_inplace(u);

    solve_lambda_section(v, u);

    delta.solve_inplace(v);

    solve_gamma_section(u, v);

    return bx;
}

DSpan1D Matrix_Corner_Block::solve_transpose_inplace(DSpan1D const bx) const
{
    assert(int(bx.extent(0)) == n);
    DSpan1D const u(bx.data_handle(), nb);
    DSpan1D const v(bx.data_handle() + nb, k);

    solve_gamma_section_transpose(v, u);

    delta.solve_transpose_inplace(v);

    solve_lambda_section_transpose(u, v);

    q_block->solve_transpose_inplace(u);

    return bx;
}

DSpan2D Matrix_Corner_Block::solve_multiple_inplace(DSpan2D const bx) const
{
    assert(int(bx.extent(0)) == n);
    for (std::size_t i(0); i < bx.extent(0); ++i) {
        DSpan1D const b(bx.data_handle() + n * i, n);
        solve_inplace(b);
    }
    return bx;
}
