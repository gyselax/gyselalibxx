// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>
#include <vector>

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

#include "ddc_aliases.hpp"

enum class BCond { PERIODIC, DIRICHLET };

/**
 * @brief A class which implements Lagrange polynomials.
 *
 * A simple class which provides the possibility to evaluate
 * an interpolation of a function known only over a restricted set of nodes.
 */
template <class Execspace, class GridInterp, BCond BcMin, BCond BcMax>
class Lagrange
{
    static_assert((BcMin == BCond::PERIODIC) == (BcMax == BCond::PERIODIC));
    using IdxInterp = Idx<GridInterp>;
    using IdxStepInterp = IdxStep<GridInterp>;
    using IdxRangeInterp = IdxRange<GridInterp>;

    using CoordDimI = Coord<typename GridInterp::continuous_dimension_type>;

private:
    IdxRangeInterp m_idx_range;
    IdxRangeInterp m_inner_idx_range;
    Field<double, IdxRangeInterp, std::experimental::layout_right, typename Execspace::memory_space>
            m_lagrange_coeffs;
    CoordDimI m_left_bound;
    CoordDimI m_right_bound;
    IdxStepInterp m_poly_support;

public:
    /**
     * @brief Usual Constructor
     *
     * @param[in] degree integer which correspond to the degree of interpolation.
     * @param[in] x_nodes_fnodes Chunkspan of nodes and associated values of the function.
     * @param[in] idx_range along interest direction, usedful in periodic case
     * @param[in] ghost DiscretVector which gives the number of ghosted points
     */
    KOKKOS_FUNCTION Lagrange(
            int degree,
            Field<double,
                  IdxRangeInterp,
                  std::experimental::layout_right,
                  typename Execspace::memory_space> x_nodes_fnodes,
            IdxRangeInterp idx_range,
            IdxStepInterp ghost)
        : m_idx_range(idx_range)
        , m_inner_idx_range(idx_range.remove(ghost, ghost))
        , m_lagrange_coeffs(x_nodes_fnodes)
        , m_left_bound(ddc::coordinate(m_inner_idx_range.front()))
        , m_right_bound(ddc::coordinate(m_inner_idx_range.back()))
        , m_poly_support(degree + 1)
    {
    }

    /**
     * @brief Evaluates the approximated value of a function on a point 
     * current values at a known set of interpolation points.
     *
     * @param[in] x_interp a node where we want to evaluate the function.
     *
     * @return The evaluation of Lagrange interpolation at the point x_intercept.
     */
    KOKKOS_FUNCTION double evaluate(CoordDimI x_interp) const;

private:
    IdxInterp getclosest(CoordDimI value) const;

    KOKKOS_FUNCTION IdxInterp getclosest_binsearch(CoordDimI value) const;

    KOKKOS_FUNCTION double evaluate_lagrange(CoordDimI x_interp) const;

    KOKKOS_FUNCTION double apply_bc(CoordDimI x_interp) const;

    KOKKOS_FUNCTION double compute_basis(
            CoordDimI x_interp,
            IdxInterp j,
            IdxRangeInterp polynom_subidx_range) const;
};

template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
Idx<GridInterp> Lagrange<Execspace, GridInterp, BcMin, BcMax>::getclosest(CoordDimI value) const
{
    assert(value >= m_left_bound && value <= m_right_bound);
    auto it = std::
            min_element(m_idx_range.begin(), m_idx_range.end(), [=](IdxInterp x, IdxInterp y) {
                return std::abs(ddc::coordinate(x) - value) < std::abs(ddc::coordinate(y) - value);
            });
    return *it;
}
template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION Idx<GridInterp> Lagrange<Execspace, GridInterp, BcMin, BcMax>::
        getclosest_binsearch(CoordDimI x_interp) const
{
    assert(x_interp >= m_left_bound && x_interp <= m_right_bound);
    IdxInterp begin = m_idx_range.front();
    IdxInterp end = m_idx_range.back();
    IdxInterp elm_cell = begin + (end - begin) / 2;
    while (x_interp < ddc::coordinate(elm_cell)
           || x_interp > ddc::coordinate(elm_cell + IdxStepInterp(1))) {
        if (x_interp < ddc::coordinate(elm_cell)) {
            end = elm_cell;
        } else {
            begin = elm_cell;
        }

        elm_cell = begin + (end - begin) / 2;
    }
    return elm_cell;
}

/**
 * @brief Computes the basis function of Lagrange.
 *
 * @param[in] x_interp a node where we want to evaluate the function.
 * @param[in] j index of the basis.
 * @param[in] polynom_subindex range a part of the mesh centered around x_interp
 *
 * @return The value of the basis at x_intercept
 */
template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, GridInterp, BcMin, BcMax>::compute_basis(
        CoordDimI x_interp,
        IdxInterp j,
        IdxRangeInterp polynom_subidx_range) const
{
    double w(1);
    CoordDimI coord_j = ddc::coordinate(j);
    for (IdxInterp const i : polynom_subidx_range) {
        CoordDimI coord_i = ddc::coordinate(i);
        if (coord_i != coord_j) {
            w *= (x_interp - coord_i) / (coord_j - coord_i);
        }
    }
    return w;
}

template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, GridInterp, BcMin, BcMax>::apply_bc(
        CoordDimI x_interp) const
{
    CoordDimI bc_val = x_interp;
    const double d = m_right_bound - m_left_bound;
    if constexpr (BcMin == BCond::PERIODIC) {
        bc_val -= std::floor((x_interp - m_left_bound) / d) * d;
    } else {
        if (x_interp < m_left_bound && BcMin == BCond::DIRICHLET) {
            bc_val = m_left_bound;
        } else if (x_interp > m_right_bound && BcMax == BCond::DIRICHLET) {
            bc_val = m_right_bound;
        }
    }
    return evaluate_lagrange(bc_val);
}

template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, GridInterp, BcMin, BcMax>::evaluate(
        CoordDimI x_interp) const
{
    if (x_interp < m_left_bound || m_right_bound < x_interp) {
        return apply_bc(x_interp);
    } else {
        return evaluate_lagrange(x_interp);
    }
}

template <typename Execspace, class GridInterp, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, GridInterp, BcMin, BcMax>::evaluate_lagrange(
        CoordDimI x_intern) const
{
    assert(x_intern >= m_left_bound && x_intern <= m_right_bound);

    IdxInterp begin, end;
    IdxInterp icell = getclosest_binsearch(x_intern);
    IdxInterp mid = icell;
    if (mid >= m_inner_idx_range.back() && BcMax == BCond::PERIODIC) {
        begin = mid - m_poly_support / 2;
        end = std::min(m_inner_idx_range.back(), begin + m_poly_support);
    } else if (mid <= m_inner_idx_range.front() && BcMin == BCond::PERIODIC) {
        begin = std::max(m_idx_range.front(), mid - m_poly_support / 2);
        end = begin + m_poly_support;
    } else {
        if (m_inner_idx_range.front() + m_poly_support / 2 > mid) {
            begin = m_inner_idx_range.front();
        } else {
            begin = mid - m_poly_support / 2;
        }
        end = std::min(m_idx_range.back() + IdxStepInterp(1), begin + m_poly_support);
    }

    if (end == m_idx_range.back())
        begin = end - m_poly_support;
    IdxStepInterp npoints_subidx_range(end - begin);
    IdxRangeInterp subidx_range(begin, npoints_subidx_range);
    double p = 0;
    for (IdxInterp const ix : subidx_range) {
        p += compute_basis(x_intern, ix, subidx_range) * m_lagrange_coeffs(ix);
    }
    return p;
}
