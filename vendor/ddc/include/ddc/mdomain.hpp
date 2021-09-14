#pragma once

#include "ddc/mcoord.hpp"
#include "ddc/mesh.hpp"
#include "ddc/non_uniform_mesh.hpp"
#include "ddc/product_mesh.hpp"
#include "ddc/rcoord.hpp"
#include "ddc/uniform_mesh.hpp"

template <class Mesh>
class MDomain
{
    static_assert(is_mesh_v<Mesh>, "A mesh is required");

public:
    using mesh_type = Mesh;

    using rcoord_type = RCoord<typename Mesh::rdim_type>;

    using mcoord_type = MCoord<Mesh>;

    struct Iterator;

private:
    /// The mesh on which this domain is defined
    mesh_type const& m_mesh;

    mcoord_type m_lbound;

    mcoord_type m_ubound;

public:
    MDomain() = default;

    template <class UboundType>
    inline constexpr MDomain(Mesh const& mesh, UboundType&& ubound) noexcept
        : m_mesh(mesh)
        , m_lbound(0)
        , m_ubound(std::forward<UboundType>(ubound))
    {
    }

    template <class LBoundType, class UboundType>
    inline constexpr MDomain(Mesh const& mesh, LBoundType&& lbound, UboundType&& ubound) noexcept
        : m_mesh(mesh)
        , m_lbound(std::forward<LBoundType>(lbound))
        , m_ubound(std::forward<UboundType>(ubound))
    {
    }

    MDomain(MDomain const& x) = default;

    MDomain(MDomain&& x) = default;

    ~MDomain() = default;

    MDomain& operator=(MDomain const& x) = default;

    MDomain& operator=(MDomain&& x) = default;

    friend constexpr bool operator==(MDomain const& xx, MDomain const& yy)
    {
        return (&xx == &yy)
               || (xx.mesh() == yy.mesh() && xx.m_lbound == yy.m_lbound
                   && xx.m_ubound == yy.m_ubound);
    }

    template <class OMesh>
    friend constexpr bool operator==(MDomain const& xx, const MDomain<OMesh>& yy)
    {
        return false;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    friend constexpr bool operator!=(MDomain const& xx, MDomain const& yy)
    {
        return !operator==(xx, yy);
    }

    template <class OMesh>
    friend constexpr bool operator!=(MDomain const& xx, const MDomain<OMesh>& yy)
    {
        return false;
    }
#endif

    inline constexpr mesh_type const& mesh() const noexcept
    {
        return m_mesh;
    }

    inline constexpr mesh_type& mesh() noexcept
    {
        return m_mesh;
    }

    inline constexpr mcoord_type front() const noexcept
    {
        return m_lbound;
    }

    inline constexpr mcoord_type back() const noexcept
    {
        return m_ubound;
    }

    inline constexpr std::size_t size() const noexcept
    {
        return m_ubound + 1 - m_lbound;
    }

    inline constexpr rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_mesh.to_real(icoord);
    }

    inline constexpr rcoord_type rmin() const noexcept
    {
        return mesh().to_real(front());
    }

    inline constexpr rcoord_type rmax() const noexcept
    {
        return mesh().to_real(back());
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
        return Iterator {front()};
    }

    constexpr Iterator cbegin() const noexcept
    {
        return begin();
    }

    constexpr Iterator end() const noexcept
    {
        return Iterator {back() + 1};
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

template <class Mesh>
struct MDomain<Mesh>::Iterator
{
private:
    typename Mesh::mcoord_type m_value = typename Mesh::mcoord_type();

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = typename Mesh::mcoord_type;

    using difference_type = MLengthElement;

    Iterator() = default;

    constexpr explicit Iterator(typename Mesh::mcoord_type __value) : m_value(__value) {}

    constexpr typename Mesh::mcoord_type operator*() const noexcept
    {
        return m_value;
    }

    constexpr Iterator& operator++()
    {
        ++m_value;
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
        --m_value;
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
            m_value += static_cast<MCoordElement>(__n);
        else
            m_value -= static_cast<MCoordElement>(-__n);
        return *this;
    }

    constexpr Iterator& operator-=(difference_type __n)
    {
        if (__n >= difference_type(0))
            m_value -= static_cast<MCoordElement>(__n);
        else
            m_value += static_cast<MCoordElement>(-__n);
        return *this;
    }

    constexpr MCoordElement operator[](difference_type __n) const
    {
        return MCoordElement(m_value + __n);
    }

    friend constexpr bool operator==(Iterator const& xx, Iterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    friend constexpr bool operator!=(Iterator const& xx, Iterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    friend constexpr bool operator<(Iterator const& xx, Iterator const& yy)
    {
        return xx.m_value < yy.m_value;
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
        return (yy.m_value > xx.m_value) ? (-static_cast<difference_type>(yy.m_value - xx.m_value))
                                         : (xx.m_value - yy.m_value);
    }
};
