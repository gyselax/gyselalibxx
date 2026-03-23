// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"

template <class BSplines>
using knot_grid_t = std::conditional_t<
        ddc::is_uniform_bsplines_v<BSplines>,
        ddc::UniformBsplinesKnots<BSplines>,
        ddc::NonUniformBsplinesKnots<BSplines>>;

template <class QuadratureGrid, class BSplines>
class QuadratureBetweenBreakPoints
{
public:
    using KnotGrid = knot_grid_t<BSplines>;

private:
    IdxRange<QuadratureGrid> m_idx_range_quadrature;
    int m_npoints_per_cell;

public:
    explicit QuadratureBetweenBreakPoints(IdxRange<QuadratureGrid> idx_range_quadrature)
        : m_idx_range_quadrature(idx_range_quadrature)
        , m_npoints_per_cell(idx_range_quadrature.size() / ddc::discrete_space<BSplines>().ncells())
    {
    }

    KOKKOS_FUNCTION IdxRange<QuadratureGrid> get_relevant_idx_range(
            IdxRange<KnotGrid> idx_range_knots) const
    {
        IdxStep<KnotGrid> cell_offset(
                idx_range_knots.front()
                - ddc::discrete_space<BSplines>().break_point_domain().front());

        if constexpr (BSplines::is_periodic) {
            if (cell_offset < 0)
                cell_offset += ddc::discrete_space<BSplines>().nbasis();
        }

        Idx<QuadratureGrid> idx_range_start
                = m_idx_range_quadrature.front() + cell_offset.value() * m_npoints_per_cell;
        IdxStep<QuadratureGrid> idx_range_result_len(
                idx_range_knots.extents().value() * m_npoints_per_cell);
        IdxRange<QuadratureGrid> idx_range_result(idx_range_start, idx_range_result_len);

        assert(idx_range_result.size() > 0);

        return idx_range_result;
    }
};

template <class Grid1D, class BSplines>
class QuadratureOnInterpolationPoints
{
public:
    using KnotGrid = knot_grid_t<BSplines>;
    using QuadratureGrid = Grid1D;

private:
    IdxRange<QuadratureGrid> m_idx_range_quadrature;
    int m_npoints_per_cell;

public:
    explicit QuadratureOnInterpolationPoints(IdxRange<Grid1D> idx_range_quadrature)
        : m_idx_range_quadrature(idx_range_quadrature)
    {
    }

    KOKKOS_FUNCTION IdxRange<QuadratureGrid> get_relevant_idx_range(
            IdxRange<KnotGrid> idx_range_knots) const
    {
        IdxStep<KnotGrid> first_cell_offset(
                idx_range_knots.front()
                - ddc::discrete_space<BSplines>().break_point_domain().front());
        IdxStep<KnotGrid> last_cell_offset(
                idx_range_knots.back() + 1
                - ddc::discrete_space<BSplines>().break_point_domain().front());

        if constexpr (BSplines::is_periodic) {
            if (first_cell_offset < 0)
                first_cell_offset += ddc::discrete_space<BSplines>().nbasis();
            if (last_cell_offset < 0)
                last_cell_offset += ddc::discrete_space<BSplines>().nbasis();
        }

        // Approximately 1 interpolation point per cell
        Idx<QuadratureGrid> idx_likely_first(
                m_idx_range_quadrature.front() + first_cell_offset.value());
        while (ddc::coordinate(idx_likely_first) >= ddc::coordinate(idx_range_knots.front())
               and idx_likely_first > m_idx_range_quadrature.front())
            [[unlikely]]
            {
                idx_likely_first -= 1;
            }
        while (ddc::coordinate(idx_likely_first) < ddc::coordinate(idx_range_knots.front())) {
            idx_likely_first += 1;
        }
        Idx<QuadratureGrid> idx_likely_last(
                m_idx_range_quadrature.front() + last_cell_offset.value()
                + ddc::discrete_space<BSplines>().degree());
        while (ddc::coordinate(idx_likely_last) > ddc::coordinate(idx_range_knots)) {
            idx_likely_last -= 1;
        }

        IdxRange<QuadratureGrid>
                idx_range_result(idx_likely_first, idx_likely_last - idx_likely_first);

        assert(idx_range_result.size() > 0);

        return idx_range_result;
    }
};
