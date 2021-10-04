#pragma once

#include <cstdint>
#include <tuple>

#include "ddc/mcoord.hpp"
#include "ddc/mdomain.hpp"
#include "ddc/mesh.hpp"
#include "ddc/rcoord.hpp"
#include "ddc/taggedtuple.hpp"

template <class... Meshes>
class ProductMDomain;

template <class... Meshes>
class ProductMDomain
{
    template <class Mesh>
    using rdim_t = typename Mesh::rdim_type;

    template <class Mesh>
    using storage_t = detail::MDomain<Mesh>;

    // static_assert((... && is_mesh_v<Meshes>), "A template parameter is not a mesh");

    static_assert((... && (Meshes::rank() == 1)), "Only rank 1 meshes are allowed.");

    template <class...>
    friend class ProductMDomain;

public:
    using rcoord_type = RCoord<rdim_t<Meshes>...>;

    using mcoord_type = MCoord<Meshes...>;

    using mlength_type = MLength<Meshes...>;

private:
    std::tuple<storage_t<Meshes>...> m_domains;

    explicit constexpr ProductMDomain(storage_t<Meshes> const&... domains) : m_domains(domains...)
    {
    }

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
        : m_domains(std::get<storage_t<Meshes>>(domains.m_domains)...)
    {
    }

    /** Construct a ProductMDomain starting from (0, ..., 0) with size points.
     * @param meshes the discrete dimensions on which the domain is constructed
     * @param size the number of points in each direction
     */
    constexpr ProductMDomain(Meshes const&... meshes, mlength_type const& size)
        : m_domains(storage_t<Meshes>(::get<Meshes>(meshes), 0, ::get<Meshes>(size))...)
    {
    }

    /** Construct a ProductMDomain starting from (0, ..., 0) with size points.
     * @param mesh the discrete space on which the domain is constructed
     * @param size the number of points in each direction
     */
    constexpr ProductMDomain(ProductMesh<Meshes...> const& mesh, mlength_type const& size)
        : m_domains(storage_t<Meshes>(::get<Meshes>(mesh), 0, ::get<Meshes>(size))...)
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
        : m_domains(storage_t<Meshes>(meshes, ::get<Meshes>(lbound), ::get<Meshes>(size))...)
    {
    }

    /** Construct a ProductMDomain starting from lbound with size points.
     * @param mesh the discrete space on which the domain is constructed
     * @param lbound the lower bound in each direction
     * @param size the number of points in each direction
     */
    constexpr ProductMDomain(
            ProductMesh<Meshes...> const& mesh,
            mcoord_type const& lbound,
            mlength_type const& size)
        : m_domains(storage_t<
                    Meshes>(::get<Meshes>(mesh), ::get<Meshes>(lbound), ::get<Meshes>(size))...)
    {
    }

    ProductMDomain(ProductMDomain const& x) = default;

    ProductMDomain(ProductMDomain&& x) = default;

    ~ProductMDomain() = default;

    ProductMDomain& operator=(ProductMDomain const& x) = default;

    ProductMDomain& operator=(ProductMDomain&& x) = default;

    constexpr bool operator==(ProductMDomain const& other) const
    {
        return m_domains == other.m_domains;
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    constexpr bool operator!=(ProductMDomain const& other) const
    {
        return !(*this == other);
    }
#endif

    template <class... QueryMeshes>
    ProductMesh<Meshes...> mesh() const
    {
        return select<QueryMeshes...>(mesh());
    }

    ProductMesh<Meshes...> mesh() const
    {
        return ProductMesh<Meshes...>(std::get<storage_t<Meshes>>(m_domains).mesh()...);
    }

    std::size_t size() const
    {
        return (1ul * ... * std::get<storage_t<Meshes>>(m_domains).size());
    }

    constexpr mlength_type extents() const noexcept
    {
        return mcoord_type(std::get<storage_t<Meshes>>(m_domains).size()...);
    }

    constexpr mcoord_type front() const noexcept
    {
        return mcoord_type(std::get<storage_t<Meshes>>(m_domains).front()...);
    }

    constexpr mcoord_type back() const noexcept
    {
        return mcoord_type(std::get<storage_t<Meshes>>(m_domains).back()...);
    }

    rcoord_type to_real(mcoord_type const& icoord) const noexcept
    {
        return rcoord_type(
                std::get<storage_t<Meshes>>(m_domains).to_real(::get<Meshes>(icoord))...);
    }

    rcoord_type rmin() const noexcept
    {
        return rcoord_type(std::get<storage_t<Meshes>>(m_domains).rmin()...);
    }

    rcoord_type rmax() const noexcept
    {
        return rcoord_type(std::get<storage_t<Meshes>>(m_domains).rmax()...);
    }

    template <class... OMeshes>
    constexpr auto restrict(ProductMDomain<OMeshes...> const& odomain) const
    {
        assert(((std::get<storage_t<OMeshes>>(m_domains).front()
                 <= std::get<storage_t<OMeshes>>(odomain.m_domains).front())
                && ...));
        assert(((std::get<storage_t<OMeshes>>(m_domains).back()
                 >= std::get<storage_t<OMeshes>>(odomain.m_domains).back())
                && ...));
        return ProductMDomain(get_slicer_for<Meshes>(odomain)...);
    }

    constexpr bool empty() const noexcept
    {
        return size() == 0;
    }

    constexpr explicit operator bool()
    {
        return !empty();
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    auto begin() const
    {
        return std::get<I>(m_domains).begin();
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    auto end() const
    {
        return std::get<I>(m_domains).end();
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    auto cbegin() const
    {
        return std::get<I>(m_domains).begin();
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    auto cend() const
    {
        return std::get<I>(m_domains).end();
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    constexpr decltype(auto) operator[](std::size_t __n)
    {
        return begin()[__n];
    }

    template <std::size_t N = sizeof...(Meshes), std::enable_if_t<N == 1, std::size_t> I = 0>
    constexpr decltype(auto) operator[](std::size_t __n) const
    {
        return begin()[__n];
    }

private:
    template <class QueryMesh, class... OMeshes>
    auto get_slicer_for(ProductMDomain<OMeshes...> const& c) const
    {
        if constexpr (in_tags_v<QueryMesh, detail::TypeSeq<OMeshes...>>) {
            return std::get<storage_t<QueryMesh>>(c.m_domains);
        } else {
            return std::get<storage_t<QueryMesh>>(m_domains);
        }
    }
};

template <class... QueryMeshes, class... Meshes>
constexpr auto select(ProductMDomain<Meshes...> const& domain)
{
    return ProductMDomain(
            select<QueryMeshes...>(domain.mesh()),
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
