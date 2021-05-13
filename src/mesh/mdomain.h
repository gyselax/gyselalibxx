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

    template <class...>
    friend class RegularMesh;

public:
    template <class... OTags>
    inline constexpr RegularMesh(const RegularMesh<OTags...>& other) noexcept
        : m_origin(other.m_origin)
        , m_step(other.m_step)
    {
    }

    template <class... OTags>
    inline constexpr RegularMesh(RegularMesh<OTags...>&& other) noexcept
        : m_origin(std::move(other.m_origin))
        , m_step(std::move(other.m_step))
    {
    }

    template <class OriginType, class StepType>
    inline constexpr RegularMesh(OriginType&& origin, StepType&& step) noexcept
        : m_origin(std::forward<OriginType>(origin))
        , m_step(std::forward<StepType>(step))
    {
    }

    friend constexpr bool operator==(const RegularMesh& xx, const RegularMesh& yy)
    {
        return (&xx == &yy) || (xx.m_origin == yy.m_origin && xx.m_step == yy.m_step);
    }

    friend constexpr bool operator!=(const RegularMesh& xx, const RegularMesh& yy)
    {
        return (&xx != &yy) && (xx.m_origin != yy.m_origin || xx.m_step != yy.m_step);
    }

    template <class... OTags>
    friend constexpr bool operator==(const RegularMesh& xx, const RegularMesh<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const RegularMesh& xx, const RegularMesh<OTags...>& yy)
    {
        return false;
    }

    static inline constexpr size_t rank() noexcept
    {
        return sizeof...(Tags);
    }

    inline constexpr RCoord_ origin() const noexcept
    {
        return m_origin;
    }

    inline constexpr RLength_ step() const noexcept
    {
        return m_step;
    }

    template <class... OTags>
    inline constexpr RCoord<OTags...> to_real(const MCoord<OTags...> icoord) const noexcept
    {
        return RCoord_((::get<OTags>(origin()) + ::get<OTags>(icoord) * ::get<OTags>(m_step))...);
    }
};

template <class... Tags>
std::ostream& operator<<(std::ostream& out, RegularMesh<Tags...> const& dom)
{
    return out << "RegularMesh( origin=" << dom.origin() << ", unitvec=" << dom.step() << " )";
}

template <class... Tags>
using Mesh = RegularMesh<Tags...>;

using MeshT = Mesh<Dim::T>;

using MeshX = Mesh<Dim::X>;

using MeshVx = Mesh<Dim::Vx>;

using MeshXVx = Mesh<Dim::X, Dim::Vx>;


namespace detail {

template <
        class... SelectedTags,
        class TagsHead,
        class... TagsQueue,
        class... AllTags,
        class SliceSpec>
inline constexpr auto append_if_all(
        RegularMesh<SelectedTags...>,
        RegularMesh<TagsHead>,
        RegularMesh<AllTags...> m,
        SliceSpec) noexcept
{
    static_assert(
            std::is_integral_v<
                    SliceSpec> || std::is_same_v<std::experimental::all_type, SliceSpec>);
    if constexpr (std::is_integral_v<SliceSpec>) {
        return RegularMesh<SelectedTags...>(m);
    } else {
        //TODO: s+h
        return RegularMesh<TagsHead, SelectedTags...>(m);
    }
}


template <class TagsHead, class... TagsQueue, class SliceSpecsHead, class... SliceSpecsQueue>
inline constexpr auto select_tags(
        RegularMesh<TagsHead, TagsQueue...> m,
        SliceSpecsHead&& h,
        SliceSpecsQueue&&... q) noexcept
{
    static_assert(sizeof...(TagsQueue) == sizeof...(SliceSpecsQueue));
    if constexpr (sizeof...(TagsQueue) > 0) {
        return append_if_all(
                select_tags(RegularMesh<TagsQueue...>(m), q...),
                RegularMesh<TagsHead>(m),
                m,
                h);
    } else {
        return append_if_all(RegularMesh<>(RCoord<>(), RLength<>()), m, m, h);
    }
}

} // namespace detail

template <class... Tags, class... SliceSpecs>
inline constexpr auto submesh(const RegularMesh<Tags...>& mesh, SliceSpecs... slices) noexcept
{
    using ReturnType = decltype(detail::select_tags(mesh, std::forward<SliceSpecs>(slices)...));
    return ReturnType(mesh);
}


template <class... Tags>
class RegularMDomain : public RegularMesh<Tags...>
{
public:
    using RegularMesh_ = RegularMesh<Tags...>;

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
