#pragma once

#include "mcoord.h"
#include "rcoord.h"
#include "uniformmesh.h"


template <class... Tags>
class RegularMDomain : public UniformMesh<Tags...>
{
public:
    using RegularMesh_ = UniformMesh<Tags...>;

    using RCoord_ = RCoord<Tags...>;

    using MCoord_ = MCoord<Tags...>;

    using Mesh = RegularMesh_;

    template <class...>
    friend class RegularMDomain;

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
    /// step size
    MCoord_ m_lbound;

    /// step size
    MCoord_ m_ubound;

public:
    template <class... OTags>
    inline constexpr RegularMDomain(const RegularMDomain<OTags...>& other) noexcept
        : RegularMesh_(other)
        , m_lbound(other.m_lbound)
        , m_ubound(other.m_ubound)
    {
    }

    template <class... OTags>
    inline constexpr RegularMDomain(RegularMDomain<OTags...>&& other) noexcept
        : RegularMesh_(std::move(other))
        , m_lbound(std::move(other.m_lbound))
        , m_ubound(std::move(other.m_ubound))
    {
    }

    template <class MeshType, class UboundType>
    inline constexpr RegularMDomain(MeshType&& mesh, UboundType&& ubound) noexcept
        : RegularMesh_(std::forward<MeshType>(mesh))
        , m_lbound(0)
        , m_ubound(std::forward<UboundType>(ubound))
    {
    }

    template <class MeshType, class LboundType, class UboundType>
    inline constexpr RegularMDomain(
            MeshType&& mesh,
            LboundType&& lbound,
            UboundType&& ubound) noexcept
        : RegularMesh_(std::forward<MeshType>(mesh))
        , m_lbound(std::forward<LboundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
        assert((m_lbound == MCoord_ {0ul}) && "non null lbound is not supported yet");
    }

    template <class OriginType, class StepType, class LboundType, class UboundType>
    inline constexpr RegularMDomain(
            OriginType&& origin,
            StepType&& step,
            LboundType&& lbound,
            UboundType&& ubound) noexcept
        : RegularMesh_(std::forward<OriginType>(origin), std::forward<StepType>(step))
        , m_lbound(std::forward<LboundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
        assert((m_lbound == MCoord_ {0ul}) && "non null lbound is not supported yet");
    }

    friend constexpr bool operator==(const RegularMDomain& xx, const RegularMDomain& yy)
    {
        return (&xx == &yy)
               || (static_cast<RegularMesh_>(xx) == static_cast<RegularMesh_>(yy)
                   && xx.m_lbound == yy.m_lbound && xx.m_ubound == yy.m_ubound);
    }

    friend constexpr bool operator!=(const RegularMDomain& xx, const RegularMDomain& yy)
    {
        return !operator==(xx, yy);
    }

    template <class... OTags>
    friend constexpr bool operator==(const RegularMDomain& xx, const RegularMDomain<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const RegularMDomain& xx, const RegularMDomain<OTags...>& yy)
    {
        return false;
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
        return this->to_real(lbound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmin() const noexcept
    {
        return this->to_real(lbound<OTags...>());
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
        return this->to_real(ubound());
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> rmax() const noexcept
    {
        return this->to_real(ubound<OTags...>());
    }

    template <class QueryTag>
    inline constexpr ptrdiff_t extent() const noexcept
    {
        return get<QueryTag>(m_ubound) - get<QueryTag>(static_cast<MCoord_>(m_lbound));
    }

    inline constexpr ptrdiff_t size() const noexcept
    {
        return ((extent<Tags>()) * ...);
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
        static_assert(sizeof...(Tags) == 1);
        return Iterator {static_cast<MCoord_>(m_lbound)};
    }

    constexpr Iterator cbegin() const noexcept
    {
        static_assert(sizeof...(Tags) == 1);
        return begin();
    }

    constexpr Iterator end() const noexcept
    {
        static_assert(sizeof...(Tags) == 1);
        return Iterator {static_cast<MCoord_>(m_ubound)};
    }

    constexpr Iterator cend() const noexcept
    {
        static_assert(sizeof...(Tags) == 1);
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
};

template <class... Tags>
std::ostream& operator<<(std::ostream& out, RegularMDomain<Tags...> const& dom)
{
    out << "RegularMDomain( origin=" << dom.origin() << ", unitvec=" << dom.step()
        << ", lower_bound=" << dom.lbound() << ", upper_bound(excluded)=" << dom.ubound() << " )";
    return out;
}

/* For now MDomain is just an alias to RegularMDomain, in the long run, we should use a tuple-based
 * solutions to have different types in each dimension
 */
template <class... Tags>
using MDomain = RegularMDomain<Tags...>;

using MDomainT = MDomain<Dim::T>;

using MDomainX = MDomain<Dim::X>;

using MDomainVx = MDomain<Dim::Vx>;

using MDomainXVx = MDomain<Dim::X, Dim::Vx>;
