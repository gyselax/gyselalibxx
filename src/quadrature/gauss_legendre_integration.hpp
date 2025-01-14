// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <vector>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"

/**
 * @brief A structure containing the weights and positions associated
 * with a Gauss-Legendre quadrature using NPoints points.
 *
 * @tparam NPoints The number of points in the quadrature scheme.
 */
template <std::size_t NPoints>
struct GaussLegendreCoefficients
{
    /**
     * @brief The positions on the domain [-1,1] where the function
     * should be evaluated.
     */
    static std::array<long double, NPoints> pos;
    /**
     * @brief The weights that the function should be multiplied by
     * in the Gauss-Legendre quadrature.
     */
    static std::array<long double, NPoints> weight;
};

template struct GaussLegendreCoefficients<1>;
template struct GaussLegendreCoefficients<2>;
template struct GaussLegendreCoefficients<3>;
template struct GaussLegendreCoefficients<4>;
template struct GaussLegendreCoefficients<5>;
template struct GaussLegendreCoefficients<6>;
template struct GaussLegendreCoefficients<7>;
template struct GaussLegendreCoefficients<8>;
template struct GaussLegendreCoefficients<9>;
template struct GaussLegendreCoefficients<10>;

/**
 * @brief An operator for constructing a Gauss-Legendre quadrature.
 * @tparam GLGrid The grid describing the Gauss-Legendre points.
 * @tparam NPoints The number of points in the Gauss-Legendre scheme
 */
template <class GLGrid, std::size_t NPoints>
class GaussLegendre
{
    static_assert(ddc::is_non_uniform_point_sampling_v<GLGrid>);

    using Dim = typename GLGrid::continuous_dimension_type;

public:
    /**
     * @brief The order of the quadrature scheme.
     * @return The order of the quadrature scheme.
     */
    static constexpr std::size_t order()
    {
        return 2 * NPoints - 1;
    }

private:
    GaussLegendreCoefficients<NPoints> m_glc;
    std::size_t const m_nbcells;
    IdxRange<GLGrid> m_valid_idx_range;
    std::vector<double> m_cell_lengths;

public:
    /**
     * @brief A constructor of the GaussLegendre class.
     * @param[in] mesh_edges_begin An iterator pointing to the first element in an
     *              iterable decscribing the edges of the cells on which the
     *              Gauss-Legendre quadrature is calculated.
     * @param[in] mesh_edges_end An iterator pointing to the end of an
     *              iterable decscribing the edges of the cells on which the
     *              Gauss-Legendre quadrature is calculated.
     */
    template <class InputIt>
    GaussLegendre(InputIt mesh_edges_begin, InputIt mesh_edges_end)
        : m_nbcells(mesh_edges_end - mesh_edges_begin - 1)
        , m_valid_idx_range(
                  Idx<GLGrid>(0),
                  IdxStep<GLGrid>((mesh_edges_end - mesh_edges_begin - 1) * NPoints))
        , m_cell_lengths(m_nbcells)
    {
        ddc::init_discrete_space<GLGrid>(get_sampling(mesh_edges_begin, mesh_edges_end));
    }

    /**
     * @brief A constructor of the GaussLegendre class.
     * @param[in] mesh_edges An initialiser list containing the edges of the cells on
     *          which the Gauss-Legendre quadrature is calculated.
     */
    explicit GaussLegendre(std::initializer_list<Coord<Dim>> const mesh_edges)
        : GaussLegendre(mesh_edges.begin(), mesh_edges.end())
    {
    }

    /**
     * @brief A constructor of the GaussLegendre class.
     * @param[in] mesh_edges A constant Field containing the edges of the cells on
     *          which the Gauss-Legendre quadrature is calculated.
     */
    template <class Grid1D>
    explicit GaussLegendre(ConstField<Coord<Dim>, IdxRange<Grid1D>, Kokkos::HostSpace> mesh_edges)
        : m_nbcells(mesh_edges.size() - 1)
        , m_valid_idx_range(Idx<GLGrid>(0), IdxStep<GLGrid>((mesh_edges.size() - 1) * NPoints))
        , m_cell_lengths(m_nbcells)
    {
        ddc::init_discrete_space<GLGrid>(get_sampling(mesh_edges));
    }

    /**
     * @brief Get the index range of the points of the Gauss-Legendre quadrature.
     * @return The index range where functions should be evaluated.
     */
    IdxRange<GLGrid> get_idx_range() const
    {
        return m_valid_idx_range;
    }

    /**
     * @brief Get a FieldMem containing the coefficients for the Gauss-Legendre quadrature.
     * @return The Gauss-Legendre quadrature.
     */
    template <class ExecSpace>
    DFieldMem<IdxRange<GLGrid>, typename ExecSpace::memory_space> gauss_legendre_coefficients()
            const
    {
        DFieldMem<IdxRange<GLGrid>, typename ExecSpace::memory_space> coefficients_alloc(
                m_valid_idx_range);
        auto coefficients_host = ddc::create_mirror_view(get_field(coefficients_alloc));
        ddc::for_each(m_valid_idx_range, [&](Idx<GLGrid> ix) {
            int i = (ix - m_valid_idx_range.front()) / NPoints;
            int j = (ix - m_valid_idx_range.front()) % NPoints;
            coefficients_host(ix) = m_cell_lengths[i] * m_glc.weight[j];
        });
        ddc::parallel_deepcopy(get_field(coefficients_alloc), coefficients_host);
        return coefficients_alloc;
    }

private:
    void get_sampling_on_cell(
            std::vector<Coord<Dim>>& grid,
            Coord<Dim> x0,
            Coord<Dim> x1,
            int i,
            int& k)
    {
        double const l = 0.5 * (ddc::get<Dim>(x1) - ddc::get<Dim>(x0));
        m_cell_lengths[i] = l;
        Coord<Dim> const c(0.5 * (x0 + x1));
        for (int j(0); j < NPoints; ++k, ++j) {
            grid[k] = c + double(l * m_glc.pos[j]);
        }
    }

    template <class InputIt>
    std::vector<Coord<Dim>> get_sampling(InputIt mesh_edges_begin, InputIt mesh_edges_end)
    {
        std::vector<Coord<Dim>> grid(m_nbcells * NPoints);

        Coord<Dim> x0;
        Coord<Dim> x1(*mesh_edges_begin);

        Idx<GLGrid> start = m_valid_idx_range.front();
        int k(0);
        for (int i(0); i < m_nbcells; ++i) {
            x0 = x1;
            x1 = *(++mesh_edges_begin);
            get_sampling_on_cell(grid, x0, x1, i, k);
        }
        return grid;
    }

    template <class Grid1D>
    std::vector<Coord<Dim>> get_sampling(
            ConstField<Coord<Dim>, IdxRange<Grid1D>, Kokkos::HostSpace> mesh_edges)
    {
        std::vector<Coord<Dim>> grid(m_nbcells * NPoints);

        IdxRange<Grid1D> idx_range = ::get_idx_range(mesh_edges);
        Idx<Grid1D> mesh_idx(idx_range.front());

        Coord<Dim> x0;
        Coord<Dim> x1 = mesh_edges(mesh_idx);

        Idx<GLGrid> start = m_valid_idx_range.front();
        int k(0);
        for (int i(0); i < m_nbcells; ++i) {
            x0 = x1;
            x1 = mesh_edges(++mesh_idx);
            get_sampling_on_cell(grid, x0, x1, i, k);
        }
        return grid;
    }
};
