#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include <Kokkos_Core.hpp>

#include "sll/matrix.hpp"
#include "sll/matrix_dense.hpp"
#include "sll/matrix_periodic_banded.hpp"

using std::max;

Matrix_Periodic_Banded::Matrix_Periodic_Banded(
        int const n,
        int const kl,
        int const ku,
        std::unique_ptr<Matrix> q)
    : Matrix_Corner_Block(n, max(kl, ku), std::move(q), max(kl, ku) + 1, max(kl, ku))
    , kl(kl)
    , ku(ku)
{
}

double Matrix_Periodic_Banded::get_element(int const i, int j) const
{
    assert(i >= 0);
    assert(i < n);
    assert(j >= 0);
    assert(i < n);
    if (i >= nb && j < nb) {
        int d = j - i;
        if (d > n / 2)
            d -= n;
        if (d < -n / 2)
            d += n;

        if (d < -kl || d > ku)
            return 0.0;
        if (d > 0)
            return lambda(j, i - nb);
        else
            return lambda(j - nb + k + 1, i - nb);
    } else {
        return Matrix_Corner_Block::get_element(i, j);
    }
}

void Matrix_Periodic_Banded::set_element(int const i, int j, double const a_ij)
{
    assert(i >= 0);
    assert(i < n);
    assert(j >= 0);
    assert(i < n);
    if (i >= nb && j < nb) {
        int d = j - i;
        if (d > n / 2)
            d -= n;
        if (d < -n / 2)
            d += n;

        if (d < -kl || d > ku) {
            assert(std::fabs(a_ij) < 1e-20);
            return;
        }

        if (d > 0)
            lambda(j, i - nb) = a_ij;
        else
            lambda(j - nb + k + 1, i - nb) = a_ij;
    } else {
        Matrix_Corner_Block::set_element(i, j, a_ij);
    }
}

void Matrix_Periodic_Banded::calculate_delta_to_factorize()
{
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            double val = 0.0;
            // Upper diagonals in lambda, lower diagonals in Abm_1_gamma
            for (int l = 0; l < i + 1; ++l) {
                val += lambda(l, i) * Abm_1_gamma(j, l);
            }
            // Lower diagonals in lambda, upper diagonals in Abm_1_gamma
            for (int l = i + 1; l < k + 1; ++l) {
                int l_full = nb - 1 - k + l;
                val += lambda(l, i) * Abm_1_gamma(j, l_full);
            }
            delta.set_element(i, j, delta.get_element(i, j) - val);
        }
    }
}

DSpan1D Matrix_Periodic_Banded::solve_lambda_section(DSpan1D const v, DView1D const u) const
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in lambda
        for (int j = 0; j <= i; ++j) {
            v(i) -= lambda(j, i) * u(j);
        }
        // Lower diagonals in lambda
        for (int j = i + 1; j < k + 1; ++j) {
            v(i) -= lambda(j, i) * u(nb - 1 - k + j);
        }
    }
    return v;
}

DSpan1D Matrix_Periodic_Banded::solve_lambda_section_transpose(DSpan1D const u, DView1D const v)
        const
{
    for (int i = 0; i < k; ++i) {
        // Upper diagonals in lambda
        for (int j = 0; j <= i; ++j) {
            u(j) -= lambda(j, i) * v(i);
        }
        // Lower diagonals in lambda
        for (int j = i + 1; j < k + 1; ++j) {
            u(nb - 1 - k + j) -= lambda(j, i) * v(i);
        }
    }
    return u;
}
