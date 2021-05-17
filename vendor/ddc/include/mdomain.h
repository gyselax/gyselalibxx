#pragma once

#include "mcoord.h"
#include "rcoord.h"
#include "uniformmesh.h"


template <class Mesh>
class MDomainImpl
{
public:
    using Mesh_ = Mesh;

    using RCoord_ = typename Mesh::RCoord_;

    using MCoord_ = typename Mesh::MCoord_;

    template <class>
    friend class MDomainImpl;

    struct Iterator
    {
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

        friend constexpr bool operator==(const Iterator& xx, const Iterator& yy)
        {
            return xx._M_value == yy._M_value;
        }

        friend constexpr bool operator!=(const Iterator& xx, const Iterator& yy)
        {
            return xx._M_value != yy._M_value;
        }

        friend constexpr bool operator<(const Iterator& xx, const Iterator& yy)
        {
            return xx._M_value < yy._M_value;
        }

        friend constexpr bool operator>(const Iterator& xx, const Iterator& yy)
        {
            return yy < xx;
        }

        friend constexpr bool operator<=(const Iterator& xx, const Iterator& yy)
        {
            return !(yy < xx);
        }

        friend constexpr bool operator>=(const Iterator& xx, const Iterator& yy)
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

        friend constexpr difference_type operator-(const Iterator& xx, const Iterator& yy)
        {
            return (yy._M_value > xx._M_value)
                           ? (-static_cast<difference_type>(yy._M_value - xx._M_value))
                           : (xx._M_value - yy._M_value);
        }

    private:
        MCoordElement _M_value = MCoordElement();
    };

private:
    /// The mesh on which this domain is defined
    Mesh_ m_mesh;

    /// step size
    MCoord_ m_lbound;

    /// step size
    MCoord_ m_ubound;

public:
    template <class OMesh>
    inline constexpr MDomainImpl(const MDomainImpl<OMesh>& other) noexcept
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
        assert((m_lbound == MCoord_ {0ul}) && "non null lbound is not supported yet");
    }

    template <class OriginType, class StepType, class LboundType, class UboundType>
    inline constexpr MDomainImpl(
            OriginType&& origin,
            StepType&& step,
            LboundType&& lbound,
            UboundType&& ubound) noexcept
        : m_mesh(std::forward<OriginType>(origin), std::forward<StepType>(step))
        , m_lbound(std::forward<LboundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
        assert((m_lbound == MCoord_ {0ul}) && "non null lbound is not supported yet");
    }

    friend constexpr bool operator==(const MDomainImpl& xx, const MDomainImpl& yy)
    {
        return (&xx == &yy)
               || (xx.mesh() == yy.mesh() && xx.m_lbound == yy.m_lbound
                   && xx.m_ubound == yy.m_ubound);
    }

    friend constexpr bool operator!=(const MDomainImpl& xx, const MDomainImpl& yy)
    {
        return !operator==(xx, yy);
    }

    template <class... OTags>
    friend constexpr bool operator==(const MDomainImpl& xx, const MDomainImpl<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const MDomainImpl& xx, const MDomainImpl<OTags...>& yy)
    {
        return false;
    }

    inline constexpr Mesh_ const& mesh() const noexcept
    {
        return m_mesh;
    }

    inline constexpr Mesh_& mesh() noexcept
    {
        return m_mesh;
    }

    inline constexpr MCoord_& lbound() noexcept
    {
        return m_lbound;
    }

    inline constexpr const MCoord_& lbound() const noexcept
    {
        return m_lbound;
    }

    template <class... OTags>
    inline constexpr MCoord<OTags...> lbound() const noexcept
    {
        return lbound();
    }

    inline constexpr RCoord_ rmin() const noexcept
    {
        return mesh().to_real(lbound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmin() const noexcept
    {
        return mesh().to_real(lbound<OTags...>());
    }

    inline constexpr MCoord_& ubound() noexcept
    {
        return m_ubound;
    }

    inline constexpr const MCoord_& ubound() const noexcept
    {
        return m_ubound;
    }

    template <class... OTags>
    inline constexpr MCoord<OTags...> ubound() const noexcept
    {
        return ubound();
    }

    inline constexpr RCoord_ rmax() const noexcept
    {
        return mesh().to_real(ubound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmax() const noexcept
    {
        return mesh().to_real(ubound<OTags...>());
    }

    template <class QueryTag>
    inline constexpr ptrdiff_t extent() const noexcept
    {
        return get<QueryTag>(m_ubound) - get<QueryTag>(static_cast<MCoord_>(m_lbound));
    }

    inline constexpr ptrdiff_t size() const noexcept
    {
        return size(std::make_index_sequence<Mesh_::rank()>());
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
        static_assert(Mesh_::rank() == 1);
        return Iterator {static_cast<MCoord_>(m_lbound)};
    }

    constexpr Iterator cbegin() const noexcept
    {
        static_assert(Mesh_::rank() == 1);
        return begin();
    }

    constexpr Iterator end() const noexcept
    {
        static_assert(Mesh_::rank() == 1);
        return Iterator {static_cast<MCoord_>(m_ubound)};
    }

    constexpr Iterator cend() const noexcept
    {
        static_assert(Mesh_::rank() == 1);
        return end();
    }

    constexpr decltype(auto) back()
    {
        assert(!empty());
        return *(--end());
    }

    constexpr decltype(auto) operator[](ptrdiff_t __n)
    {
        return begin()[__n];
    }

    constexpr decltype(auto) operator[](ptrdiff_t __n) const
    {
        return begin()[__n];
    }

private:
    template <size_t... Idxs>
    inline constexpr ptrdiff_t size(std::index_sequence<Idxs...>) const noexcept
    {
        return ((m_ubound[Idxs] - m_lbound[Idxs]) * ...);
    }
};

template <class Mesh>
std::ostream& operator<<(std::ostream& out, MDomainImpl<Mesh> const& dom)
{
    out << "MDomain( origin=" << dom.mesh().origin() << ", unitvec=" << dom.mesh().step()
        << ", lower_bound=" << dom.lbound() << ", upper_bound(excluded)=" << dom.ubound() << " )";
    return out;
}

/* For now MDomain is just an alias to MDomain, in the long run, we should use a tuple-based
 * solutions to have different types in each dimension
 */
template <class... Tags>
using MDomain = MDomainImpl<UniformMesh<Tags...>>;

template <class... Tags>
using RegularMDomain = MDomainImpl<UniformMesh<Tags...>>;

template <class... Tags>
using UniformMDomain = MDomainImpl<UniformMesh<Tags...>>;

using MDomainX = UniformMDomain<Dim::X>;

using MDomainVx = UniformMDomain<Dim::Vx>;

using MDomainXVx = UniformMDomain<Dim::X, Dim::Vx>;
