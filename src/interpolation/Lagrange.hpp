// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>
#include <vector>

#include <ddc/ddc.hpp>


enum class BCond { PERIODIC, DIRICHLET };

/**
 * @brief A class which implements Lagrange polynomials.
 *
 * A simple class which provides the possibility to evaluate
 * an interpolation of a function known only over a restricted set of nodes.
 */
template <class DDim, BCond BcMin, BCond BcMax>
class Lagrange
{
    static_assert((BcMin == BCond::PERIODIC) == (BcMax == BCond::PERIODIC));

    ddc::DiscreteDomain<DDim> m_domain;
    int m_deg;
    int m_n;
    int m_ghost;
    ddc::Chunk<double, ddc::DiscreteDomain<DDim>> m_Chunk;
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
    Lagrange(
            int order,
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> x_nodes_fnodes,
            ddc::DiscreteDomain<DDim> domain,
            ddc::DiscreteVector<DDim> ghost)
        : m_domain(domain)
        , m_deg(order)
        , m_n(x_nodes_fnodes.domain().size())
        , m_ghost(ghost.value())
        , m_Chunk(domain)
    {
        if constexpr (BcMin == BCond::PERIODIC) {
            for (int k = 0; k < m_ghost; k++) {
                m_Chunk(m_domain.front() + ddc::DiscreteVector<DDim>(k))
                        = x_nodes_fnodes(ddc::DiscreteElement<DDim>(m_n + k));
            }

            ddc::for_each(
                    ddc::policies::serial_host,
                    x_nodes_fnodes.domain(),
                    [&](ddc::DiscreteElement<DDim> const ix) { m_Chunk(ix) = x_nodes_fnodes(ix); });
            assert(domain.back() == ddc::DiscreteElement<DDim>(m_n - 1 + 2 * m_ghost));

            for (int k = 0; k <= m_ghost; k++) {
                m_Chunk(ddc::DiscreteElement<DDim>(m_n - 1 + m_ghost + k))
                        = x_nodes_fnodes(ddc::DiscreteElement<DDim>(m_ghost + k));
            }
        } else {
            ddc::deepcopy(m_Chunk, x_nodes_fnodes);
        }
        m_left_bound = ddc::coordinate(ddc::DiscreteElement<DDim>(m_ghost));
        m_right_bound = ddc::coordinate(ddc::DiscreteElement<DDim>(m_n + m_ghost - 1));
    };

    /**
     * @brief Evaluates the approximated value of a function on a point 
     * current values at a known set of interpolation points.
     *
     * @param[in] x_interp a node where we want to evaluate the function.
     *
     * @return The evaluation of Lagrange interpolation at the point x_intercept.
     */
    double evaluate(double x_interp);

    /**
     * @brief Computes the basis function of Lagrange.
     *
     * @param[in] x_interp a node where we want to evaluate the function.
     * @param[in] j index of the basis.
     * @param[in] polynom_subdomain a part of the mesh centered around x_interp
     *
     * @return The value of the basis at x_intercept
     */
    double compute_basis(
            double x_interp,
            ddc::DiscreteElement<DDim> j,
            ddc::DiscreteDomain<DDim> polynom_subdomain);

private:
    auto getclosest(double value) const;

    int getclosest_binsearch(double value) const;

    double evaluate_lagrange(double x_interp);

    double apply_bc(double x_interp);
};

template <class DDim, BCond BcMin, BCond BcMax>
auto Lagrange<DDim, BcMin, BcMax>::getclosest(double value) const
{
    assert(value >= m_left_bound && value <= m_right_bound);
    auto it = std::min_element(m_domain.begin(), m_domain.end(), [=](auto x, auto y) {
        return std::abs(ddc::coordinate(x) - value) < std::abs(ddc::coordinate(y) - value);
    });
    return *it - m_Chunk.domain().front();
}

template <class DDim, BCond BcMin, BCond BcMax>
int Lagrange<DDim, BcMin, BcMax>::getclosest_binsearch(double x_interp) const
{
    int begin = m_domain.front().uid();
    int end = m_domain.back().uid();
    int icell = (begin + end) / 2;

    while (x_interp < ddc::coordinate(ddc::DiscreteElement<DDim>(icell))
           || x_interp > ddc::coordinate(ddc::DiscreteElement<DDim>(icell + 1))) {
        if (x_interp < ddc::coordinate((ddc::DiscreteElement<DDim>(icell)))) {
            end = icell;
        } else {
            begin = icell;
        }

        icell = (begin + end) / 2;
    }

    return icell;
}

template <class DDim, BCond BcMin, BCond BcMax>
double Lagrange<DDim, BcMin, BcMax>::compute_basis(
        double x_interp,
        ddc::DiscreteElement<DDim> j,
        ddc::DiscreteDomain<DDim> polynom_subdomain)
{
    return ddc::transform_reduce(
            ddc::policies::serial_host,
            polynom_subdomain,
            1.0,
            ddc::reducer::prod<double>(),
            [&](ddc::DiscreteElement<DDim> const ix) {
                ddc::Coordinate<typename DDim::continuous_dimension_type> coord_i
                        = ddc::coordinate(ix);
                ddc::Coordinate<typename DDim::continuous_dimension_type> coord_j
                        = ddc::coordinate(j);
                if (coord_i != coord_j) {
                    return (x_interp - coord_i) / (coord_j - coord_i);
                } else {
                    return 1.0;
                }
            });
}

template <class DDim, BCond BcMin, BCond BcMax>
double Lagrange<DDim, BcMin, BcMax>::apply_bc(double x_interp)
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

template <class DDim, BCond BcMin, BCond BcMax>
double Lagrange<DDim, BcMin, BcMax>::evaluate(double x_interp)
{
    double interpol_val = x_interp;
    if (x_interp < m_left_bound || m_right_bound < x_interp) {
        interpol_val = apply_bc(x_interp);
    } else {
        interpol_val = evaluate_lagrange(x_interp);
    }

    return interpol_val;
}

template <class DDim, BCond BcMin, BCond BcMax>
double Lagrange<DDim, BcMin, BcMax>::evaluate_lagrange(double x_intern)
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
    ddc::DiscreteElement<DDim> lbound_subdom(begin);
    ddc::DiscreteVector<DDim> npoints_subdomain(end - begin);
    ddc::DiscreteDomain<DDim> subdomain(lbound_subdom, npoints_subdomain);

    double p = 0;
    ddc::for_each(ddc::policies::serial_host, subdomain, [&](ddc::DiscreteElement<DDim> const ix) {
        p += compute_basis(x_intern, ix, subdomain) * m_Chunk(ix);
    });

    return p;
}
