#pragma once

#include "mcoord.h"
#include "rcoord.h"

template <class... Tags>
class RegularMesh
{
public:
    using RCoord_ = RCoord<Tags...>;

    using RLength_ = RLength<Tags...>;

    using MCoord_ = MCoord<Tags...>;

private:
    /// origin
    RCoord_ m_origin;

    /// step size
    RLength_ m_step;

public:
    inline constexpr RegularMesh ( RCoord_ origin, RLength_ step ) noexcept
        : m_origin ( origin )
        , m_step ( step )
    {
    }

    inline constexpr RCoord_ origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr RLength_ step() const noexcept
    {
        return m_step;
    }

    inline constexpr RCoord_ rcoord ( const MCoord_ icoord ) const noexcept
    {
        return {origin() + icoord * m_step};
    }

    inline constexpr RCoord_ operator() ( const MCoord_ icoord ) const noexcept
    {
        return rcoord ( icoord );
    }

    static inline constexpr size_t rank () noexcept
    {
        return sizeof... ( Tags );
    }
};


template <class... Tags>
class RegularMDomain: public RegularMesh<Tags...>
{
public:
    using RegularMesh_ = RegularMesh<Tags...>;

    using RCoord_ = RCoord<Tags...>;

    using MCoord_ = MCoord<Tags...>;

    using Mesh = RegularMesh_;

private:
    /// step size
    MCoord_ m_begin;

    /// step size
    MCoord_ m_end;

public:
    inline constexpr RegularMDomain ( RegularMesh_ mesh, RCoord_ begin, RCoord_ end ) noexcept
        : RegularMesh_ ( std::move ( mesh ) )
        , m_begin ( std::move ( begin ) )
        , m_end ( std::move ( end ) )
    {
    }

    inline constexpr RegularMDomain ( RegularMesh_ mesh, RCoord_ end ) noexcept
        : RegularMesh_ ( std::move ( mesh ) )
        , m_begin ( 0 )
        , m_end ( std::move ( end ) )
    {
    }

    inline constexpr RegularMDomain ( RCoord_ origin, double step, RCoord_ begin, RCoord_ end ) noexcept
        : RegularMesh_ ( std::move ( origin ), step )
        , m_begin ( std::move ( begin ) )
        , m_end ( std::move ( end ) )
    {
    }

    inline constexpr RegularMDomain ( RCoord_ origin, double step, RCoord_ end ) noexcept
        : RegularMesh_ ( std::move ( origin ), step )
        , m_begin ( 0 )
        , m_end ( std::move ( end ) )
    {
    }

    inline constexpr MCoord_& begin() noexcept
    {
        return m_begin;
    }

    inline constexpr const MCoord_& begin() const noexcept
    {
        return m_begin;
    }

    inline constexpr const MCoord_& cbegin() const noexcept
    {
        return m_begin;
    }

    inline constexpr MCoord_& end() noexcept
    {
        return m_end;
    }

    inline constexpr const MCoord_& end() const noexcept
    {
        return m_end;
    }

    inline constexpr const MCoord_& cend() const noexcept
    {
        return m_end;
    }
};


/* For now MDomain is just an alias to RegularMDomain, in the long run, we should use a tuple-based
 * solutions to have different types in each dimension
 */
template <class... Tags>
using MDomain = RegularMDomain<Tags...>;

using MDomainX = MDomain<Dim::X>;

using MDomainVx = MDomain<Dim::Vx>;

using MDomainXVx = MDomain<Dim::X, Dim::Vx>;
