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

private:
    ddc::DiscreteDomain<DDimI> m_domain;
    int m_deg;
    int m_n;
    int m_ghost;
    ddc::ChunkSpan<
            double,
            ddc::DiscreteDomain<DDimI>,
            std::experimental::layout_right,
            typename Execspace::memory_space>
            m_ChunkSpan;
    double m_left_bound;
    double m_right_bound;

public:
    /**
     * @brief Usual Constructor
     *
     * @param[in] order integer which correspond to the degree of interpolation.
     * @param[in] x_nodes_fnodes Chunkspan of nodes and associated values of the function.
     * @param[in] domain along interest direction, usedful in periodic case
     * @param[in] ghost DiscretVector which gives the number of ghosted points
     */
    KOKKOS_FUNCTION Lagrange(
            int order,
            ddc::ChunkSpan<
                    double,
                    ddc::DiscreteDomain<DDimI>,
                    std::experimental::layout_right,
                    typename Execspace::memory_space> x_nodes_fnodes,
            ddc::DiscreteDomain<DDimI> domain,
            ddc::DiscreteVector<DDimI> ghost)
        : m_domain(domain)
        , m_deg(order)
        , m_n(x_nodes_fnodes.domain().size())
        , m_ghost(ghost.value())
        , m_ChunkSpan(x_nodes_fnodes)
    {
        m_left_bound = ddc::coordinate(ddc::DiscreteElement<DDimI>(m_ghost));
        m_right_bound = ddc::coordinate(ddc::DiscreteElement<DDimI>(m_n + m_ghost - 1));
    };

    /**
     * @brief Evaluates the approximated value of a function on a point 
     * current values at a known set of interpolation points.
     *
     * @param[in] x_interp a node where we want to evaluate the function.
     *
     * @return The evaluation of Lagrange interpolation at the point x_intercept.
     */
    KOKKOS_FUNCTION double evaluate(double x_interp) const;

private:
    auto getclosest(double value) const;

    KOKKOS_FUNCTION int getclosest_binsearch(double value) const;

    KOKKOS_FUNCTION double evaluate_lagrange(double x_interp) const;

    KOKKOS_FUNCTION double apply_bc(double x_interp) const;

    KOKKOS_FUNCTION double compute_basis(
            double x_interp,
            ddc::DiscreteElement<DDimI> j,
            ddc::DiscreteDomain<DDimI> polynom_subdomain) const;
};

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
auto Lagrange<Execspace, DDimI, BcMin, BcMax>::getclosest(double value) const
{
    assert(value >= m_left_bound && value <= m_right_bound);
    auto it = std::min_element(m_domain.begin(), m_domain.end(), [=](auto x, auto y) {
        return std::abs(ddc::coordinate(x) - value) < std::abs(ddc::coordinate(y) - value);
    });
    return *it - m_ChunkSpan.domain().front();
}
template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION int Lagrange<Execspace, DDimI, BcMin, BcMax>::getclosest_binsearch(
        double x_interp) const
{
    std::size_t begin = m_domain.front().uid();
    std::size_t end = m_domain.back().uid();
    std::size_t icell = (begin + end) / 2;

    while (x_interp < ddc::coordinate(ddc::DiscreteElement<DDimI>(icell))
           || x_interp > ddc::coordinate(ddc::DiscreteElement<DDimI>(icell + 1))) {
        if (x_interp < ddc::coordinate((ddc::DiscreteElement<DDimI>(icell)))) {
            end = icell;
        } else {
            begin = icell;
        }

        icell = (begin + end) / 2;
    }

    return icell;
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
        double x_interp,
        ddc::DiscreteElement<DDimI> j,
        ddc::DiscreteDomain<DDimI> polynom_subdomain) const
{
    double w(1);
    ddc::Coordinate<typename DDimI::continuous_dimension_type> coord_j = ddc::coordinate(j);
    for (ddc::DiscreteElement<DDimI> const i : polynom_subdomain) {
        ddc::Coordinate<typename DDimI::continuous_dimension_type> coord_i = ddc::coordinate(i);
        if (coord_i != coord_j) {
            w *= (x_interp - coord_i) / (coord_j - coord_i);
        }
    }
    return w;
}

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::apply_bc(
        double x_interp) const
{
    double bc_val = x_interp;
    const double d = std::abs(m_right_bound - m_left_bound);
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
        double x_interp) const
{
    double interpol_val = x_interp;
    if (x_interp < m_left_bound || m_right_bound < x_interp) {
        interpol_val = apply_bc(x_interp);
    } else {
        interpol_val = evaluate_lagrange(x_interp);
    }

    return interpol_val;
}

template <typename Execspace, class DDimI, BCond BcMin, BCond BcMax>
KOKKOS_INLINE_FUNCTION double Lagrange<Execspace, DDimI, BcMin, BcMax>::evaluate_lagrange(
        double x_intern) const
{
    int begin, end;
    int mid = getclosest_binsearch(x_intern);

    if (mid >= int(m_domain.size()) - 1 - 2 * m_ghost && BcMax == BCond::PERIODIC) {
        begin = mid - int(0.5 * m_deg);
        end = std::min(int(m_domain.size()) - 1, begin + m_deg);
    } else if (mid <= m_ghost + 1 && BcMin == BCond::PERIODIC) {
        begin = std::max(0, mid - int(0.5 * m_deg));
        end = begin + m_deg;
    } else {
        begin = std::max(m_ghost, mid - int(0.5 * m_deg));
        end = std::min(m_n, begin + m_deg + 1);
    }

    if (end == int(m_domain.size()) - 1)
        begin = end - m_deg;

    ddc::DiscreteElement<DDimI> lbound_subdom(begin);
    ddc::DiscreteVector<DDimI> npoints_subdomain(end - begin);
    ddc::DiscreteDomain<DDimI> subdomain(lbound_subdom, npoints_subdomain);
    double p = 0;
    for (ddc::DiscreteElement<DDimI> const ix : subdomain) {
        p += compute_basis(x_intern, ix, subdomain) * m_ChunkSpan(ix);
    }

    return p;
}
