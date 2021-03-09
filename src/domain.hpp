#pragma once

#include <array>
#include <tuple>

template <class T, class D>
class Size {
    T m_value;

public:
    inline Size(T value)
        : m_value(value)
    {
    }

    inline Size& operator=(T value)
    {
        m_value = value;
        return *this;
    }

    inline operator T() const { return m_value; }
};

template <class T, class D>
class Coordinate {
public:
    using Diff = T;

private:
    Diff m_value;

public:
    inline Coordinate(Diff value)
        : m_value(value)
    {
    }

    inline Coordinate& operator=(Diff value)
    {
        m_value = value;
        return *this;
    }

    inline operator Diff() const { return m_value; }
};

template <class C>
/// A semi-open interval [lbound, ubound)
class Interval {
public:
    using Coordinate = C;

    using Diff = typename Coordinate::Diff;

private:
    std::array<Coordinate, 2> m_bounds;

public:
    inline Interval(std::array<Coordinate, 2> value)
        : m_bounds(value)
    {
    }

    inline Interval(Coordinate lbound, Coordinate ubound)
        : m_bounds { lbound, ubound }
    {
    }

    inline Interval& operator=(std::array<Coordinate, 2> value)
    {
        m_bounds = value;
        return *this;
    }

    inline operator std::array<Coordinate, 2>() const { return m_bounds; }

    inline Coordinate lbound() const { return m_bounds[0]; }

    inline Interval& lbound(Coordinate lb)
    {
        m_bounds[0] = lb;
        return *this;
    }

    inline Coordinate ubound() const { return m_bounds[1]; }

    inline Interval& ubound(Coordinate ub)
    {
        m_bounds[1] = ub;
        return *this;
    }

    inline Diff size() const { return m_bounds[1]; }

    inline Diff& size() { return m_bounds[1]; }
};

template <class... D>
class RealDomain {
public:
    using Diff = std::tuple<Size<double, D>...>;

    using GlobalCoordinate = std::tuple<Coordinate<Size<double, D>, D>...>;

    using LocalCoordinate = std::tuple<Coordinate<Size<double, D>, D>...>;

    using GlobalInterval = std::tuple<typename Coordinate<Size<double, D>, D>::Interval...>;

    using LocalInterval = std::tuple<typename Coordinate<Size<double, D>, D>::Interval...>;

private:
    static std::array<bool, sizeof...(D)> s_periodic;

    static Diff s_size;

    Diff m_size;

public:
};

template <class... D>
class DiscreteDomain {
public:
    using RealDomain = RealDomain<D...>;

    using Diff = std::tuple<Size<int, D>...>;

    using GlobalCoordinate = std::tuple<Coordinate<Size<int, D>, D>...>;

    using LocalCoordinate = std::tuple<Coordinate<Size<int, D>, D>...>;

    using GlobalInterval = std::tuple<typename Coordinate<Size<int, D>, D>::Interval...>;

    using LocalInterval = std::tuple<typename Coordinate<Size<int, D>, D>::Interval...>;

private:
    static std::array<double, sizeof...(D)> s_step;

public:
};

template <class RD>
struct RealDomainToDiscrete;

template <class... D>
struct RealDomainToDiscrete<RealDomain<D...>> {
    using type = DiscreteDomain<D...>;
};

class DX {
};

class DVX {
};

using RealPhaseDomain = RealDomain<DX, DVX>;

using RealSpaceDomain = RealDomain<DX>;

// using DiscretePhaseDomain = DiscreteDomain<DX, DVX>;
using DiscretePhaseDomain = RealDomainToDiscrete<RealPhaseDomain>;

using DiscreteSpaceDomain = DiscreteDomain<DX>;
