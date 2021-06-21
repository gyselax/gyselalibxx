#pragma once

#include "mcoord.h"
#include "nonuniformmesh.h"
#include "rcoord.h"
#include "uniformmesh.h"
#include "view.h"

template <class Mesh>
class MDomainImpl
{
public:
    using mesh_type = Mesh;

    using rcoord_type = typename Mesh::rcoord_type;

    using mcoord_type = typename Mesh::mcoord_type;

    template <class>
    friend class MDomainImpl;

    struct Iterator
    {
    private:
        MCoordElement _M_value = MCoordElement();

    public:
        using iterator_category = std::random_access_iterator_tag;

        using value_type = MCoordElement;

        using difference_type = MLengthElement;

        Iterator() = default;

        constexpr explicit Iterator(MCoordElement __value) : _M_value(__value) {}

        constexpr MCoordElement operator*() const noexcept
        {
            return _M_value;
        }

        constexpr Iterator& operator++()
        {
            ++_M_value;
            return *this;
        }

        constexpr Iterator operator++(int)
        {
            auto __tmp = *this;
            ++*this;
            return __tmp;
        }

        constexpr Iterator& operator--()
        {
            --_M_value;
            return *this;
        }

        constexpr Iterator operator--(int)
        {
            auto __tmp = *this;
            --*this;
            return __tmp;
        }

        constexpr Iterator& operator+=(difference_type __n)
        {
            if (__n >= difference_type(0))
                _M_value += static_cast<MCoordElement>(__n);
            else
                _M_value -= static_cast<MCoordElement>(-__n);
            return *this;
        }

        constexpr Iterator& operator-=(difference_type __n)
        {
            if (__n >= difference_type(0))
                _M_value -= static_cast<MCoordElement>(__n);
            else
                _M_value += static_cast<MCoordElement>(-__n);
            return *this;
        }

        constexpr MCoordElement operator[](difference_type __n) const
        {
            return MCoordElement(_M_value + __n);
        }

        friend constexpr bool operator==(Iterator const& xx, Iterator const& yy)
        {
            return xx._M_value == yy._M_value;
        }

        friend constexpr bool operator!=(Iterator const& xx, Iterator const& yy)
        {
            return xx._M_value != yy._M_value;
        }

        friend constexpr bool operator<(Iterator const& xx, Iterator const& yy)
        {
            return xx._M_value < yy._M_value;
        }

        friend constexpr bool operator>(Iterator const& xx, Iterator const& yy)
        {
            return yy < xx;
        }

        friend constexpr bool operator<=(Iterator const& xx, Iterator const& yy)
        {
            return !(yy < xx);
        }

        friend constexpr bool operator>=(Iterator const& xx, Iterator const& yy)
        {
            return !(xx < yy);
        }

        friend constexpr Iterator operator+(Iterator __i, difference_type __n)
        {
            return __i += __n;
        }

        friend constexpr Iterator operator+(difference_type __n, Iterator __i)
        {
            return __i += __n;
        }

        friend constexpr Iterator operator-(Iterator __i, difference_type __n)
        {
            return __i -= __n;
        }

        friend constexpr difference_type operator-(Iterator const& xx, Iterator const& yy)
        {
            return (yy._M_value > xx._M_value)
                           ? (-static_cast<difference_type>(yy._M_value - xx._M_value))
                           : (xx._M_value - yy._M_value);
        }
    };

private:
    /// The mesh on which this domain is defined
    mesh_type m_mesh;

    mcoord_type m_lbound;

    mcoord_type m_ubound;

public:
    template <class OMesh>
    inline constexpr MDomainImpl(MDomainImpl<OMesh> const& other) noexcept
        : m_mesh(other.m_mesh)
        , m_lbound(other.m_lbound)
        , m_ubound(other.m_ubound)
    {
    }

    template <class OMesh>
    inline constexpr MDomainImpl(MDomainImpl<OMesh>&& other) noexcept
        : m_mesh(std::move(other.m_mesh))
        , m_lbound(std::move(other.m_lbound))
        , m_ubound(std::move(other.m_ubound))
    {
    }

    template <class MeshType, class UboundType>
    inline constexpr MDomainImpl(MeshType&& mesh, UboundType&& ubound) noexcept
        : m_mesh(std::forward<MeshType>(mesh))
        , m_lbound(0)
        , m_ubound(std::forward<UboundType>(ubound))
    {
    }

    template <class MeshType, class LboundType, class UboundType>
    inline constexpr MDomainImpl(MeshType&& mesh, LboundType&& lbound, UboundType&& ubound) noexcept
        : m_mesh(std::forward<MeshType>(mesh))
        , m_lbound(std::forward<LboundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
        //         cannot assert in constexpr :/
        //         assert((m_lbound == mcoord_type {0ul}) && "non null lbound is not supported yet");
    }

    template <class OriginType, class StepType, class LboundType, class UboundType>
    inline constexpr MDomainImpl(
            OriginType&& rmin,
            StepType&& rmax,
            LboundType&& lbound,
            UboundType&& ubound) noexcept
        : m_mesh(
                rmin + lbound * (rmin - rmax) / (ubound - lbound),
                ((rmax - rmin) / (ubound - lbound)))
        , m_lbound(std::forward<LboundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
        //         cannot assert in constexpr :/
        //         assert((m_lbound == mcoord_type {0ul}) && "non null lbound is not supported yet");
    }

    friend constexpr bool operator==(MDomainImpl const& xx, MDomainImpl const& yy)
    {
        return (&xx == &yy)
               || (xx.mesh() == yy.mesh() && xx.m_lbound == yy.m_lbound
                   && xx.m_ubound == yy.m_ubound);
    }

    friend constexpr bool operator!=(MDomainImpl const& xx, MDomainImpl const& yy)
    {
        return !operator==(xx, yy);
    }

    template <class OMesh>
    friend constexpr bool operator==(MDomainImpl const& xx, const MDomainImpl<OMesh>& yy)
    {
        return false;
    }

    template <class OMesh>
    friend constexpr bool operator!=(MDomainImpl const& xx, const MDomainImpl<OMesh>& yy)
    {
        return false;
    }

    inline constexpr mesh_type const& mesh() const noexcept
    {
        return m_mesh;
    }

    inline constexpr mesh_type& mesh() noexcept
    {
        return m_mesh;
    }

    inline constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_mesh.to_real(icoord);
    }

    inline constexpr mcoord_type& lbound() noexcept
    {
        return m_lbound;
    }

    inline constexpr const mcoord_type& lbound() const noexcept
    {
        return m_lbound;
    }

    template <class... OTags>
    inline constexpr MCoord<OTags...> lbound() const noexcept
    {
        return lbound();
    }

    inline constexpr rcoord_type rmin() const noexcept
    {
        return mesh().to_real(lbound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmin() const noexcept
    {
        return mesh().to_real(lbound<OTags...>());
    }

    inline constexpr mcoord_type& ubound() noexcept
    {
        return m_ubound;
    }

    inline constexpr const mcoord_type& ubound() const noexcept
    {
        return m_ubound;
    }

    template <class... OTags>
    inline constexpr MCoord<OTags...> ubound() const noexcept
    {
        return ubound();
    }

    inline constexpr rcoord_type rmax() const noexcept
    {
        return mesh().to_real(ubound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmax() const noexcept
    {
        return mesh().to_real(ubound<OTags...>());
    }

    template <class QueryTag>
    inline constexpr std::size_t extent() const noexcept
    {
        return get<QueryTag>(m_ubound) - get<QueryTag>(static_cast<mcoord_type>(m_lbound));
    }

    inline constexpr ExtentsND<Mesh::rank()> extents() const noexcept
    {
        return extents(std::make_index_sequence<Mesh::rank()>());
    }

    inline constexpr std::size_t size() const noexcept
    {
        return size(std::make_index_sequence<mesh_type::rank()>());
    }

    inline constexpr bool empty() const noexcept
    {
        return size() == 0;
    }

    constexpr explicit operator bool()
    {
        return !empty();
    }

    constexpr Iterator begin() const noexcept
    {
        static_assert(mesh_type::rank() == 1);
        return Iterator {static_cast<mcoord_type>(m_lbound)};
    }

    constexpr Iterator cbegin() const noexcept
    {
        static_assert(mesh_type::rank() == 1);
        return begin();
    }

    constexpr Iterator end() const noexcept
    {
        static_assert(mesh_type::rank() == 1);
        return Iterator {static_cast<mcoord_type>(m_ubound)};
    }

    constexpr Iterator cend() const noexcept
    {
        static_assert(mesh_type::rank() == 1);
        return end();
    }

    constexpr decltype(auto) back()
    {
        assert(!empty());
        return *(--end());
    }

    constexpr decltype(auto) operator[](std::size_t __n)
    {
        return begin()[__n];
    }

    constexpr decltype(auto) operator[](std::size_t __n) const
    {
        return begin()[__n];
    }

private:
    template <size_t... Idxs>
    inline constexpr std::size_t size(std::index_sequence<Idxs...>) const noexcept
    {
        return ((m_ubound[Idxs] - m_lbound[Idxs]) * ...);
    }

    template <size_t... Idxs>
    inline constexpr ExtentsND<sizeof...(Idxs)> extents(std::index_sequence<Idxs...>) const noexcept
    {
        return ExtentsND<sizeof...(Idxs)>((m_ubound[Idxs] - m_lbound[Idxs])...);
    }
};

template <class QTag, class... CTags>
auto get_slicer_for(MCoord<CTags...> const& c)
{
    if constexpr (has_tag_v<QTag, MCoord<CTags...>>) {
        return c.template get<QTag>();
    } else {
        return std::experimental::all;
    }
}

template <class Mesh>
std::ostream& operator<<(std::ostream& out, MDomainImpl<Mesh> const& dom)
{
    out << "MDomain( origin=" << dom.mesh().origin() << ", unitvec=" << dom.mesh().step()
        << ", lower_bound=" << dom.lbound() << ", upper_bound(excluded)=" << dom.ubound() << " )";
    return out;
}

template <class... Tags>
using UniformMDomain = MDomainImpl<UniformMesh<Tags...>>;

template <class... Tag>
using NonUniformMDomain = MDomainImpl<NonUniformMesh<Tag...>>;
