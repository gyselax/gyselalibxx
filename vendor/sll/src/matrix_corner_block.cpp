#include <cassert>
#include <utility>

#include <experimental/mdspan>

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
    memset(lambda_ptr.get(), 0, sizeof(double) * k * nb);
    memset(Abm_1_gamma_ptr.get(), 0, sizeof(double) * k * nb);
}

Matrix_Corner_Block::Matrix_Corner_Block(
        int const n,
        int const k,
        std::unique_ptr<Matrix> q,
        int const lambda_size)
    : Matrix(n)
    , k(k)
    , nb(n - k)
    , Abm_1_gamma_ptr(std::make_unique<double[]>(k * nb))
    , lambda_ptr(std::make_unique<double[]>(lambda_size))
    , q_block(std::move(q))
    , delta(k)
    , Abm_1_gamma(Abm_1_gamma_ptr.get(), k, nb)
    , lambda(lambda_ptr.get(), nb, k)
{
    assert(n > 0);
    assert(k >= 0);
    assert(k <= n);
    memset(lambda_ptr.get(), 0, sizeof(double) * lambda_size);
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
    } else if (i >= nb) {
        return Abm_1_gamma(i - nb, j);
    } else {
        return lambda(i, j - nb);
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
    } else if (i >= nb) {
        Abm_1_gamma(i - nb, j) = a_ij;
    } else {
        lambda(i, j - nb) = a_ij;
    }
}

void Matrix_Corner_Block::calculate_delta_to_factorize()
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in lambda
        for (int l = 0; l < nb; ++l) {
            double const lambda_il = lambda(l, i);
            for (int j = 0; j < k; ++j) {
                double const new_val(delta.get_element(j, i) - lambda_il * Abm_1_gamma(j, l));
                delta.set_element(j, i, new_val);
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

DSpan2D Matrix_Corner_Block::solve_lambda_section(DSpan2D const v, DView2D const u) const
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in lambda
        for (int j = 0; j < nb; ++j) {
            for (int col = 0; col < v.extent(1); ++col) {
                v(i, col) -= lambda(j, i) * u(j, col);
            }
        }
    }
    return v;
}

void Matrix_Corner_Block::solve_lambda_section_transpose(DSpan1D const u, DSpan1D const v) const
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in (*lambda)
        for (int j = 0; j < nb; ++j) {
            u(j) -= lambda(j, i) * v(i);
        }
    }
}

DSpan1D Matrix_Corner_Block::solve_inplace(DSpan1D const bx) const
{
    assert(bx.extent(0) == n);
    DSpan1D const u(bx.data(), nb);
    DSpan1D const v(bx.data() + nb, k);

    q_block->solve_inplace(u);

    solve_lambda_section(v, u);

    delta.solve_inplace(v);

    for (int i = 0; i < nb; ++i) {
        double val = 0.;
        for (int j = 0; j < k; ++j) {
            val += Abm_1_gamma(j, i) * v(j);
        }
        u(i) -= val;
    }
    return bx;
}

DSpan1D Matrix_Corner_Block::solve_transpose_inplace(DSpan1D const bx) const
{
    assert(bx.extent(0) == n);
    DSpan1D const u(bx.data(), nb);
    DSpan1D const v(bx.data() + nb, k);

    delta.solve_inplace(v);

    solve_lambda_section_transpose(u, v);

    q_block->solve_inplace(u);

    for (int j = 0; j < k; ++j) {
        double val = 0.;
        for (int i = 0; i < nb; ++i) {
            val += Abm_1_gamma(j, i) * v(i);
        }
        v(j) -= val;
    }
    return bx;
}

DSpan2D Matrix_Corner_Block::solve_inplace_matrix(DSpan2D const bx) const
{
    assert(bx.extent(0) == n);
    DSpan2D const u(bx.data(), nb, bx.extent(1));
    DSpan2D const v(bx.data() + nb * bx.extent(1), k, bx.extent(1));

    q_block->solve_inplace_matrix(u);

    solve_lambda_section(v, u);

    delta.solve_inplace_matrix(v);

    for (int col = 0; col < bx.extent(1); ++col) {
        for (int i = 0; i < nb; ++i) {
            double val = 0.;
            for (int j = 0; j < k; ++j) {
                val += Abm_1_gamma(j, i) * v(j, col);
            }
            u(i, col) -= val;
        }
    }
    return bx;
}
