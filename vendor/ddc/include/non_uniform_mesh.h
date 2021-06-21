#pragma once

#include <vector>

#include <view.h>

#include "mcoord.h"
#include "rcoord.h"

template <class Tag>
class NonUniformMesh
{
public:
    using rcoord_type = RCoord<Tag>;

    using rlength_type = RLength<Tag>;

    using mcoord_type = MCoord<Tag>;

    using tag_type = Tag;

    // The two Mesh and Mesh_ need better names to avoid ambiguity
    using Mesh_ = NonUniformMesh<Tag>;

    template <class OTag>
    using Mesh = NonUniformMesh<OTag>;

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
    /// mesh points
    std::vector<rcoord_type> m_points;

    mcoord_type m_lbound;

public:
    inline constexpr NonUniformMesh(NonUniformMesh const& other) = default;

    inline constexpr NonUniformMesh(NonUniformMesh&& other) = default;

    inline constexpr NonUniformMesh(std::vector<rcoord_type>&& points, mcoord_type lbound)
        : m_points(std::move(points))
        , m_lbound(lbound)
    {
    }

    template <class InputIterable>
    inline constexpr NonUniformMesh(InputIterable const& points, mcoord_type lbound)
        : m_points(points.begin(), points.end())
        , m_lbound(lbound)
    {
    }

    inline constexpr NonUniformMesh(View1D<const rcoord_type> points, mcoord_type lbound)
        : m_points(points.data(), points.data() + points.extent(0))
        , m_lbound(lbound)
    {
    }

    template <class InputIt>
    inline constexpr NonUniformMesh(InputIt points_begin, InputIt points_end, mcoord_type lbound)
        : m_points(points_begin, points_end)
        , m_lbound(lbound)
    {
    }

    static inline constexpr size_t rank() noexcept
    {
        return 1;
    }

    inline constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_points[icoord];
    }

    friend constexpr bool operator==(NonUniformMesh const& xx, NonUniformMesh const& yy)
    {
        return xx.m_lbound == yy.m_lbound && xx.m_points == yy.m_points;
    }

    friend constexpr bool operator!=(NonUniformMesh const& xx, NonUniformMesh const& yy)
    {
        return !operator==(xx, yy);
    }

    template <class OTag>
    friend constexpr bool operator==(NonUniformMesh const& xx, NonUniformMesh<OTag> const& yy)
    {
        return false;
    }

    template <class OTag>
    friend constexpr bool operator!=(NonUniformMesh const& xx, NonUniformMesh<OTag> const& yy)
    {
        return false;
    }

    inline constexpr NonUniformMesh const& mesh() const noexcept
    {
        return *this;
    }

    inline constexpr mcoord_type lbound() const noexcept
    {
        return m_lbound;
    }

    inline constexpr rcoord_type rmin() const noexcept
    {
        return mesh().to_real(lbound());
    }

    inline constexpr mcoord_type ubound() const noexcept
    {
        return lbound() + m_points.size();
    }

    inline constexpr rcoord_type rmax() const noexcept
    {
        return mesh().to_real(ubound());
    }

    inline constexpr std::size_t size() const noexcept
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

    constexpr decltype(auto) operator[](std::size_t __n)
    {
        return begin()[__n];
    }

    constexpr decltype(auto) operator[](std::size_t __n) const
    {
        return begin()[__n];
    }
};

template <class Tag>
std::ostream& operator<<(std::ostream& out, NonUniformMesh<Tag> const& dom)
{
    return out << "NonUniformMesh(  )";
}
