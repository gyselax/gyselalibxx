#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <ddc/ddc.hpp>

#include "sll/view.hpp"

struct GaussLegendreCoefficients
{
    static constexpr std::size_t max_order = 10u;

    static constexpr std::size_t nb_coefficients = max_order * (max_order + 1) / 2;

    static std::array<long double, nb_coefficients> weight;

    static std::array<long double, nb_coefficients> pos;
};

template <class Dim>
class GaussLegendre
{
    using glc = GaussLegendreCoefficients;

public:
    static std::size_t max_order()
    {
        return glc::max_order;
    }

    explicit GaussLegendre(std::size_t n)
    {
        assert(n > 0);
        assert(n <= glc::max_order);

        std::size_t const offset = n * (n - 1) / 2;
        m_wx.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            m_wx[i].first = static_cast<double>(glc::weight[offset + i]);
            m_wx[i].second = static_cast<ddc::Coordinate<Dim>>(glc::pos[offset + i]);
        }
    }

    GaussLegendre(GaussLegendre const& x) = default;

    GaussLegendre(GaussLegendre&& x) = default;

    ~GaussLegendre() = default;

    GaussLegendre& operator=(GaussLegendre const& x) = default;

    GaussLegendre& operator=(GaussLegendre&& x) = default;

    std::size_t order() const
    {
        return m_wx.size();
    }

    template <class F>
    double integrate(F&& f, double x0, double x1) const
    {
        static_assert(std::is_invocable_r_v<double, F, double>, "Functor F not handled");
        assert(x0 <= x1);
        double const l = 0.5 * (x1 - x0);
        double const c = 0.5 * (x0 + x1);
        double integral = 0;
        for (std::pair<double, ddc::Coordinate<Dim>> const& wx : m_wx) {
            integral += wx.first * f(l * wx.second + c);
        }
        return l * integral;
    }

    template <class Domain>
    void compute_points(ddc::ChunkSpan<ddc::Coordinate<Dim>, Domain> points, double x0, double x1)
            const
    {
        assert(x0 <= x1);
        assert(points.size() == m_wx.size());
        // map the interval [-1,1] into the interval [a,b].
        double const l = 0.5 * (x1 - x0);
        double const c = 0.5 * (x0 + x1);
        int const dom_start = points.domain().front().uid();
        for (auto i : points.domain()) {
            points(i) = ddc::Coordinate<Dim>(l * m_wx[i.uid() - dom_start].second + c);
        }
    }

    template <class Domain>
    void compute_points_and_weights(
            ddc::ChunkSpan<ddc::Coordinate<Dim>, Domain> points,
            ddc::ChunkSpan<double, Domain> weights,
            ddc::Coordinate<Dim> const& x0,
            ddc::Coordinate<Dim> const& x1) const
    {
        assert(x0 <= x1);
        assert(points.size() == m_wx.size());
        // map the interval [-1,1] into the interval [a,b].
        double const l = 0.5 * (get<Dim>(x1) - get<Dim>(x0));
        ddc::Coordinate<Dim> const c = 0.5 * (x0 + x1);
        int const dom_start = points.domain().front().uid();
        for (auto i : points.domain()) {
            weights(i) = l * m_wx[i.uid() - dom_start].first;
            points(i) = l * m_wx[i.uid() - dom_start].second + c;
        }
    }

    template <class Domain1, class Domain2>
    void compute_points_and_weights_on_mesh(
            ddc::ChunkSpan<ddc::Coordinate<Dim>, Domain1> points,
            ddc::ChunkSpan<double, Domain1> weights,
            ddc::ChunkSpan<const ddc::Coordinate<Dim>, Domain2> mesh_edges) const
    {
        [[maybe_unused]] int const nbcells = mesh_edges.size() - 1;
        int const npts_gauss = m_wx.size();
        assert(points.size() == m_wx.size() * nbcells);
        assert(weights.size() == m_wx.size() * nbcells);

        ddc::for_each(
                mesh_edges.domain().remove_last(typename Domain2::mlength_type(1)),
                [&](auto icell) {
                    ddc::Coordinate<Dim> const x0 = mesh_edges(icell);
                    ddc::Coordinate<Dim> const x1 = mesh_edges(icell + 1);

                    Domain1 local_domain(
                            typename Domain1::discrete_element_type(icell.uid() * npts_gauss),
                            typename Domain1::mlength_type(npts_gauss));

                    compute_points_and_weights(points[local_domain], weights[local_domain], x0, x1);
                });
    }

private:
    std::vector<std::pair<double, ddc::Coordinate<Dim>>> m_wx;
};
