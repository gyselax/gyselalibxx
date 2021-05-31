#pragma once

#include <vector>

#include <view.h>

#include "mcoord.h"
#include "rcoord.h"

template <class Tag>
class NonUniformMesh
{
public:
    using RCoord_ = RCoord<Tag>;

    using RLength_ = RLength<Tag>;

    using MCoord_ = MCoord<Tag>;

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
    };

private:
    /// mesh points
    std::vector<RCoord_> m_points;

    MCoord_ m_lbound;

public:
    inline constexpr NonUniformMesh(NonUniformMesh const& other) = default;

    inline constexpr NonUniformMesh(NonUniformMesh const&& other) = default;

    inline constexpr NonUniformMesh(std::vector<RCoord_>&& points, MCoord_ lbound)
        : m_points(std::move(points))
        , m_lbound(lbound)
    {
    }

    template <class InputIterable>
    inline constexpr NonUniformMesh(InputIterable const& points, MCoord_ lbound)
        : m_points(points.begin(), points.end())
        , m_lbound(lbound)
    {
    }

    inline constexpr NonUniformMesh(View1D<const RCoord_> points, MCoord_ lbound)
        : m_points(points.data(), points.data() + points.extent(0))
        , m_lbound(lbound)
    {
    }

    template <class InputIt>
    inline constexpr NonUniformMesh(InputIt points_begin, InputIt points_end, MCoord_ lbound)
        : m_points(points_begin, points_end)
        , m_lbound(lbound)
    {
    }

    static inline constexpr size_t rank() noexcept
    {
        return 1;
    }

    inline constexpr MCoord_ extents() noexcept
    {
        return m_points.size();
    }

    inline constexpr RCoord_ to_real(MCoord_ const icoord) const noexcept
    {
        return m_points[icoord];
    }

    friend constexpr bool operator==(const NonUniformMesh& xx, const NonUniformMesh& yy)
    {
        return xx.m_lbound == yy.m_lbound && xx.m_points == yy.m_points;
    }

    friend constexpr bool operator!=(const NonUniformMesh& xx, const NonUniformMesh& yy)
    {
        return !operator==(xx, yy);
    }

    template <class... OTags>
    friend constexpr bool operator==(const NonUniformMesh& xx, const NonUniformMesh<OTags...>& yy)
    {
        return false;
    }

    template <class... OTags>
    friend constexpr bool operator!=(const NonUniformMesh& xx, const NonUniformMesh<OTags...>& yy)
    {
        return false;
    }

    inline constexpr NonUniformMesh const& mesh() const noexcept
    {
        return *this;
    }

    inline constexpr NonUniformMesh& mesh() noexcept
    {
        return *this;
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
        return lbound() + m_points.size();
    }

    inline constexpr const MCoord_& ubound() const noexcept
    {
        return lbound() + m_points.size();
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
        return m_points.size();
    }

    inline constexpr ptrdiff_t size() const noexcept
    {
        return m_points.size();
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
        return Iterator(lbound());
    }

    constexpr Iterator cbegin() const noexcept
    {
        return begin();
    }

    constexpr Iterator end() const noexcept
    {
        return Iterator(ubound());
    }

    constexpr Iterator cend() const noexcept
    {
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

template <class Tag>
std::ostream& operator<<(std::ostream& out, NonUniformMesh<Tag> const& dom)
{
    return out << "NonUniformMesh(  )";
}
