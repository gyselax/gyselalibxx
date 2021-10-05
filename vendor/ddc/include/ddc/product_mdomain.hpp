#pragma once

#include <cstdint>
#include <tuple>

#include "ddc/detail/product_mesh.hpp"
#include "ddc/mcoord.hpp"
#include "ddc/mesh.hpp"
#include "ddc/rcoord.hpp"
#include "ddc/taggedtuple.hpp"
#include "ddc/type_seq.hpp"

template <class Mesh>
struct ProductMDomainIterator;

template <class... Meshes>
class ProductMDomain;

template <class... Meshes>
class ProductMDomain
{
    template <class Mesh>
    using rdim_t = typename Mesh::rdim_type;

    // static_assert((... && is_mesh_v<Meshes>), "A template parameter is not a mesh");

    static_assert((... && (Meshes::rank() == 1)), "Only rank 1 meshes are allowed.");

    template <class...>
    friend class ProductMDomain;

public:
    using rcoord_type = RCoord<rdim_t<Meshes>...>;

    using mcoord_type = MCoord<Meshes...>;

    using mlength_type = MLength<Meshes...>;

private:
    detail::ProductMesh<Meshes...> m_mesh;

    MCoord<Meshes...> m_lbound;

    MCoord<Meshes...> m_ubound;

public:
    static constexpr std::size_t rank()
    {
        return (0 + ... + Meshes::rank());
    }

    ProductMDomain() = default;

    // Use SFINAE to disambiguate with the copy constructor.
    // Note that SFINAE may be redundant because a template constructor should not be selected as a copy constructor.
    template <std::size_t N = sizeof...(Meshes), class = std::enable_if_t<(N != 1)>>
    explicit constexpr ProductMDomain(ProductMDomain<Meshes> const&... domains)
        : m_mesh(domains.template mesh<Meshes>()...)
        , m_lbound(domains.front()...)
        , m_ubound(domains.back()...)
    {
    }

    /** Construct a ProductMDomain starting from (0, ..., 0) with size points.
     * @param meshes the discrete dimensions on which the domain is constructed
     * @param size the number of points in each direction
     */
    constexpr ProductMDomain(Meshes const&... meshes, mlength_type const& size)
        : m_mesh(meshes...)
        , m_lbound((get<Meshes>(size) - get<Meshes>(size))...) // Hack to have expansion of zero
        , m_ubound((get<Meshes>(size) - 1)...)
    {
    }

    /** Construct a ProductMDomain starting from lbound with size points.
     * @param meshes the discrete dimensions on which the domain is constructed
     * @param lbound the lower bound in each direction
     * @param size the number of points in each direction
     */
    constexpr ProductMDomain(
            Meshes const&... meshes,
            mcoord_type const& lbound,
            mlength_type const& size)
        : m_mesh(meshes...)
        , m_lbound(lbound)
        , m_ubound((get<Meshes>(lbound) + get<Meshes>(size) - 1)...)
    {
    }

    ProductMDomain(ProductMDomain const& x) = default;

    ProductMDomain(ProductMDomain&& x) = default;

    ~ProductMDomain() = default;

    ProductMDomain& operator=(ProductMDomain const& x) = default;

    ProductMDomain& operator=(ProductMDomain&& x) = default;

    constexpr bool operator==(ProductMDomain const& other) const
    {
        return (((get<Meshes>(m_mesh) == get<Meshes>(other.m_mesh)) && ...)
                && m_lbound == other.m_lbound && m_ubound == other.m_ubound);
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(ProductMDomain const& other) const
    {
        return !(*this == other);
    }
#endif

    template <class QueryMesh>
    auto const& mesh() const
    {
        return get<QueryMesh>(m_mesh);
    }

    std::size_t size() const
    {
        return (1ul * ... * (get<Meshes>(m_ubound) + 1 - get<Meshes>(m_lbound)));
    }

    constexpr mlength_type extents() const noexcept
    {
        return mlength_type((get<Meshes>(m_ubound) + 1 - get<Meshes>(m_lbound))...);
    }

    constexpr mcoord_type front() const noexcept
    {
        return m_lbound;
    }

    constexpr mcoord_type back() const noexcept
    {
        return m_ubound;
    }

    rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return m_mesh.to_real(icoord);
    }

    rcoord_type rmin() const noexcept
    {
        return to_real(front());
    }

    rcoord_type rmax() const noexcept
    {
        return to_real(back());
    }

    template <class... OMeshes>
    constexpr auto restrict(ProductMDomain<OMeshes...> const& odomain) const
    {
        assert(((get<OMeshes>(m_lbound) <= get<OMeshes>(odomain.m_lbound)) && ...));
        assert(((get<OMeshes>(m_ubound) >= get<OMeshes>(odomain.m_ubound)) && ...));
        const MCoord<Meshes...> myextents = extents();
        const MCoord<OMeshes...> oextents = odomain.extents();
        return ProductMDomain(
                get<Meshes>(m_mesh)...,
                MCoord<Meshes...>((get_or<Meshes>(odomain.m_lbound, get<Meshes>(m_lbound)))...),
                MCoord<Meshes...>((get_or<Meshes>(oextents, get<Meshes>(myextents)))...));
    }

    constexpr bool empty() const noexcept
    {
        return size() == 0;
    }

    constexpr explicit operator bool()
    {
        return !empty();
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    auto begin() const
    {
        return ProductMDomainIterator<Mesh0>(front());
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    auto end() const
    {
        return ProductMDomainIterator<Mesh0>(back() + 1);
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    auto cbegin() const
    {
        return ProductMDomainIterator<Mesh0>(front());
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    auto cend() const
    {
        return ProductMDomainIterator<Mesh0>(back() + 1);
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    constexpr decltype(auto) operator[](std::size_t __n)
    {
        return begin()[__n];
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    constexpr decltype(auto) operator[](std::size_t __n) const
    {
        return begin()[__n];
    }
};

template <class... QueryMeshes, class... Meshes>
constexpr ProductMDomain<QueryMeshes...> select(ProductMDomain<Meshes...> const& domain)
{
    return ProductMDomain<QueryMeshes...>(
            domain.template mesh<QueryMeshes>()...,
            select<QueryMeshes...>(domain.front()),
            select<QueryMeshes...>(domain.extents()));
}

template <class... QueryMeshes, class... Meshes>
constexpr MCoord<QueryMeshes...> extents(ProductMDomain<Meshes...> const& domain) noexcept
{
    return MCoord<QueryMeshes...>(select<QueryMeshes>(domain).size()...);
}

template <class... QueryMeshes, class... Meshes>
constexpr MCoord<QueryMeshes...> front(ProductMDomain<Meshes...> const& domain) noexcept
{
    return MCoord<QueryMeshes...>(select<QueryMeshes>(domain).front()...);
}

template <class... QueryMeshes, class... Meshes>
constexpr MCoord<QueryMeshes...> back(ProductMDomain<Meshes...> const& domain) noexcept
{
    return MCoord<QueryMeshes...>(select<QueryMeshes>(domain).back()...);
}

template <class... QueryMeshes, class... Meshes>
RCoord<QueryMeshes...> to_real(
        ProductMDomain<Meshes...> const& domain,
        MCoord<QueryMeshes...> const& icoord) noexcept
{
    return RCoord<QueryMeshes...>(
            select<QueryMeshes>(domain).to_real(select<QueryMeshes>(icoord))...);
}

template <class... QueryMeshes, class... Meshes>
RCoord<QueryMeshes...> rmin(ProductMDomain<Meshes...> const& domain) noexcept
{
    return RCoord<QueryMeshes...>(select<QueryMeshes>(domain).rmin()...);
}

template <class... QueryMeshes, class... Meshes>
RCoord<QueryMeshes...> rmax(ProductMDomain<Meshes...> const& domain) noexcept
{
    return RCoord<QueryMeshes...>(select<QueryMeshes>(domain).rmax()...);
}

namespace detail {

template <class QueryMeshesSeq>
struct Selection;

template <class... QueryMeshes>
struct Selection<detail::TypeSeq<QueryMeshes...>>
{
    template <class Domain>
    static constexpr auto select(Domain const& domain)
    {
        return ::select<QueryMeshes...>(domain);
    }
};

} // namespace detail

template <class QueryMeshesSeq, class... Meshes>
constexpr auto select_by_type_seq(ProductMDomain<Meshes...> const& domain)
{
    return detail::Selection<QueryMeshesSeq>::select(domain);
}

template <class Mesh>
struct ProductMDomainIterator
{
private:
    typename Mesh::mcoord_type m_value = typename Mesh::mcoord_type();

public:
    using iterator_category = std::random_access_iterator_tag;

    using value_type = typename Mesh::mcoord_type;

    using difference_type = MLengthElement;

    ProductMDomainIterator() = default;

    constexpr explicit ProductMDomainIterator(typename Mesh::mcoord_type __value) : m_value(__value)
    {
    }

    constexpr typename Mesh::mcoord_type operator*() const noexcept
    {
        return m_value;
    }

    constexpr ProductMDomainIterator& operator++()
    {
        ++m_value;
        return *this;
    }

    constexpr ProductMDomainIterator operator++(int)
    {
        auto __tmp = *this;
        ++*this;
        return __tmp;
    }

    constexpr ProductMDomainIterator& operator--()
    {
        --m_value;
        return *this;
    }

    constexpr ProductMDomainIterator operator--(int)
    {
        auto __tmp = *this;
        --*this;
        return __tmp;
    }

    constexpr ProductMDomainIterator& operator+=(difference_type __n)
    {
        if (__n >= difference_type(0))
            m_value += static_cast<MCoordElement>(__n);
        else
            m_value -= static_cast<MCoordElement>(-__n);
        return *this;
    }

    constexpr ProductMDomainIterator& operator-=(difference_type __n)
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

    friend constexpr bool operator==(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return xx.m_value == yy.m_value;
    }

    friend constexpr bool operator!=(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return xx.m_value != yy.m_value;
    }

    friend constexpr bool operator<(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return xx.m_value < yy.m_value;
    }

    friend constexpr bool operator>(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return yy < xx;
    }

    friend constexpr bool operator<=(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return !(yy < xx);
    }

    friend constexpr bool operator>=(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return !(xx < yy);
    }

    friend constexpr ProductMDomainIterator operator+(
            ProductMDomainIterator __i,
            difference_type __n)
    {
        return __i += __n;
    }

    friend constexpr ProductMDomainIterator operator+(
            difference_type __n,
            ProductMDomainIterator __i)
    {
        return __i += __n;
    }

    friend constexpr ProductMDomainIterator operator-(
            ProductMDomainIterator __i,
            difference_type __n)
    {
        return __i -= __n;
    }

    friend constexpr difference_type operator-(
            ProductMDomainIterator const& xx,
            ProductMDomainIterator const& yy)
    {
        return (yy.m_value > xx.m_value) ? (-static_cast<difference_type>(yy.m_value - xx.m_value))
                                         : (xx.m_value - yy.m_value);
    }
};
