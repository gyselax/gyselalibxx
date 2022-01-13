#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "sll/view.hpp"

struct GaussLegendreCoefficients
{
    static constexpr std::size_t max_order = 10u;

    static constexpr std::size_t nb_coefficients = max_order * (max_order + 1) / 2;

    static std::array<long double, nb_coefficients> weight;

    static std::array<long double, nb_coefficients> pos;
};

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
            m_wx[i].second = static_cast<double>(glc::pos[offset + i]);
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
        for (std::pair<double, double> const& wx : m_wx) {
            integral += wx.first * f(l * wx.second + c);
        }
        return l * integral;
    }

    void compute_points(Span1D<double> points, double x0, double x1) const
    {
        assert(x0 <= x1);
        assert(points.extent(0) == m_wx.size());
        // map the interval [-1,1] into the interval [a,b].
        double const l = 0.5 * (x1 - x0);
        double const c = 0.5 * (x0 + x1);
        for (std::size_t i = 0; i < m_wx.size(); ++i) {
            points(i) = l * m_wx[i].second + c;
        }
    }

    void compute_points_and_weights(
            Span1D<double> points,
            Span1D<double> weights,
            double x0,
            double x1) const
    {
        assert(x0 <= x1);
        assert(points.extent(0) == m_wx.size());
        // map the interval [-1,1] into the interval [a,b].
        double const l = 0.5 * (x1 - x0);
        double const c = 0.5 * (x0 + x1);
        for (std::size_t i = 0; i < m_wx.size(); ++i) {
            weights(i) = l * m_wx[i].first;
            points(i) = l * m_wx[i].second + c;
        }
    }

private:
    std::vector<std::pair<double, double>> m_wx;
};
