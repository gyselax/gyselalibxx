// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>
#include <vector>

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

enum class BCond { PERIODIC, DIRICHLET };

/**
 * @brief A class which implements Lagrange polynomials.
 *
 * A simple class which provides the possibility to evaluate
 * an interpolation of a function known only over a restricted set of nodes.
 */
template <class Execspace, class DDimI, BCond BcMin, BCond BcMax>
class Lagrange
{
    static_assert((BcMin == BCond::PERIODIC) == (BcMax == BCond::PERIODIC));
    using DElemI = ddc::DiscreteElement<DDimI>;
    using DVectI = ddc::DiscreteVector<DDimI>;

    using CoordDimI = ddc::Coordinate<typename DDimI::continuous_dimension_type>;

private:
    ddc::DiscreteDomain<DDimI> m_domain;
    ddc::DiscreteDomain<DDimI> m_inner_domain;
    ddc::ChunkSpan<
            double,
            ddc::DiscreteDomain<DDimI>,
            std::experimental::layout_right,
            typename Execspace::memory_space>
            m_ChunkSpan;
    CoordDimI m_left_bound;
    CoordDimI m_right_bound;
    DVectI m_poly_support;

public:
    /**
     * @brief Usual Constructor
     *
     * @param[in] degree integer which correspond to the degree of interpolation.
     * @param[in] x_nodes_fnodes Chunkspan of nodes and associated values of the function.
     * @param[in] domain along interest direction, usedful in periodic case
     * @param[in] ghost DiscretVector which gives the number of ghosted points
     */
    KOKKOS_FUNCTION Lagrange(
            int degree,
            ddc::ChunkSpan<
                    double,
                    ddc::DiscreteDomain<DDimI>,
                    std::experimental::layout_right,
                    typename Execspace::memory_space> x_nodes_fnodes,
            ddc::DiscreteDomain<DDimI> domain,
            DVectI ghost)
        : m_domain(domain)
        , m_inner_domain(domain.remove(ghost, ghost))
        , m_ChunkSpan(x_nodes_fnodes)
        , m_left_bound(ddc::coordinate(m_inner_domain.front()))
        , m_right_bound(ddc::coordinate(m_inner_domain.back()))
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
    DElemI getclosest(CoordDimI value) const;

    KOKKOS_FUNCTION DElemI getclosest_binsearch(CoordDimI value) const;

    KOKKOS_FUNCTION double evaluate_lagrange(CoordDimI x_interp) const;

    KOKKOS_FUNCTION double apply_bc(CoordDimI x_interp) const;

    KOKKOS_FUNCTION double compute_basis(
            CoordDimI x_interp,
            DElemI j,
            ddc::DiscreteDomain<DDimI> polynom_subdomain) const;
};

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
ddc::DiscreteElement<DDimI> Lagrange<Execspace, DDimI, BcMin, BcMax>::getclosest(
        CoordDimI value) const
{
    assert(value >= m_left_bound && value <= m_right_bound);
    auto it = std::min_element(m_domain.begin(), m_domain.end(), [=](DElemI x, DElemI y) {
        return std::abs(ddc::coordinate(x) - value) < std::abs(ddc::coordinate(y) - value);
    });
    return *it;
}
template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION ddc::DiscreteElement<DDimI> Lagrange<Execspace, DDimI, BcMin, BcMax>::
        getclosest_binsearch(CoordDimI x_interp) const
{
    assert(x_interp >= m_left_bound && x_interp <= m_right_bound);
    DElemI begin = m_domain.front();
    DElemI end = m_domain.back();
    DElemI elm_cell = begin + (end - begin) / 2;
    while (x_interp < ddc::coordinate(elm_cell)
           || x_interp > ddc::coordinate(elm_cell + DVectI(1))) {
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
 * @param[in] polynom_subdomain a part of the mesh centered around x_interp
 *
 * @return The value of the basis at x_intercept
 */
template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::compute_basis(
        CoordDimI x_interp,
        DElemI j,
        ddc::DiscreteDomain<DDimI> polynom_subdomain) const
{
    double w(1);
    CoordDimI coord_j = ddc::coordinate(j);
    for (DElemI const i : polynom_subdomain) {
        CoordDimI coord_i = ddc::coordinate(i);
        if (coord_i != coord_j) {
            w *= (x_interp - coord_i) / (coord_j - coord_i);
        }
    }
    return w;
}

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::apply_bc(
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

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::evaluate(
        CoordDimI x_interp) const
{
    if (x_interp < m_left_bound || m_right_bound < x_interp) {
        return apply_bc(x_interp);
    } else {
        return evaluate_lagrange(x_interp);
    }
}

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::evaluate_lagrange(
        CoordDimI x_intern) const
{
    assert(x_intern >= m_left_bound && x_intern <= m_right_bound);

    DElemI begin, end;
    DElemI icell = getclosest_binsearch(x_intern);
    DElemI mid = icell;
    if (mid >= m_inner_domain.back() && BcMax == BCond::PERIODIC) {
        begin = mid - m_poly_support / 2;
        end = std::min(m_inner_domain.back(), begin + m_poly_support);
    } else if (mid <= m_inner_domain.front() && BcMin == BCond::PERIODIC) {
        begin = std::max(m_domain.front(), mid - m_poly_support / 2);
        end = begin + m_poly_support;
    } else {
        if (m_inner_domain.front() + m_poly_support / 2 > mid) {
            begin = m_inner_domain.front();
        } else {
            begin = mid - m_poly_support / 2;
        }
        end = std::min(m_domain.back() + DVectI(1), begin + m_poly_support);
    }

    if (end == m_domain.back())
        begin = end - m_poly_support;
    DVectI npoints_subdomain(end - begin);
    ddc::DiscreteDomain<DDimI> subdomain(begin, npoints_subdomain);
    double p = 0;
    for (DElemI const ix : subdomain) {
        p += compute_basis(x_intern, ix, subdomain) * m_ChunkSpan(ix);
    }
    return p;
}
