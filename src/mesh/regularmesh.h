#pragma once

#include "mcoord.h"
#include "rcoord.h"

template <class Dimension>
class RegularMesh
{
public:
    using RCoord = RCoord<Dimension>;

    using MCoord = MCoord<Dimension>;

private:
    /// origin
    RCoord m_origin;

    /// step size
    double m_step;

public:
    inline constexpr RegularMesh(RCoord origin, double step) noexcept
        : m_origin(origin)
        , m_step(step)
    {
    }

    inline constexpr double origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr double step() const noexcept
    {
        return m_step;
    }

    inline constexpr RCoord operator()(const MCoord icoord) const noexcept
    {
        return {origin() + icoord[0] * m_step};
    }
};

template <class Dimension>
class RegularDomain: public RegularMesh<Dimension>
{
public:
    using RCoord = RCoord<Dimension>;

    using MCoord = MCoord<Dimension>;

private:
    /// origin
    RegularMesh m_mesh;

    /// step size
    RCoord m_begin;

    /// step size
    RCoord m_end;

public:
    inline constexpr RegularMesh(RegularMesh mesh, RCoord begin, RCoord end) noexcept
        : m_mesh(std::move(mesh))
        , m_begin(std::move(begin))
        , m_end(std::move(end))
    {
    }
    
    inline constexpr RegularMesh(RegularMesh mesh, RCoord end) noexcept
        : m_mesh(std::move(mesh))
        , m_begin(0)
        , m_end(std::move(end))
    {
    }

    inline constexpr RCoord& begin() noexcept
    {
        return m_begin;
    }

    inline constexpr const RCoord& begin() const noexcept
    {
        return m_begin;
    }

    inline constexpr const RCoord& cbegin() const noexcept
    {
        return m_begin;
    }

    inline constexpr RCoord& end() noexcept
    {
        return m_end;
    }

    inline constexpr const RCoord& end() const noexcept
    {
        return m_end;
    }

    inline constexpr const RCoord& cend() const noexcept
    {
        return m_end;
    }

    inline constexpr RCoord operator()(const MCoord icoord) const noexcept
    {
        return {origin() + icoord[0] * m_step};
    }
};
