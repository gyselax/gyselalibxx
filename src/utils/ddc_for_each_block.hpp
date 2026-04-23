#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>

#include <ddc/ddc.hpp>

namespace ddcHelper {

namespace detail {

constexpr bool is_power_of_2(std::size_t const n) noexcept
{
    return n > 0 && !(n & (n - 1));
}

void distribute_blocks(
        std::size_t nb_blocks,
        std::span<ddc::DiscreteVectorElement const> sizes,
        std::span<ddc::DiscreteVectorElement> nb_blocks_per_dim)
{
    assert(sizes.size() == nb_blocks_per_dim.size());

    if (!is_power_of_2(nb_blocks)) {
        throw std::runtime_error("DDC distribute_blocks expects a power of 2.");
    }

    for (ddc::DiscreteVectorElement& blocks : nb_blocks_per_dim) {
        blocks = 1;
    }

    std::size_t dim = 0;
    while (nb_blocks != 1) {
        if (sizes[dim] >= nb_blocks_per_dim[dim] * 2) {
            nb_blocks_per_dim[dim] *= 2;
            nb_blocks /= 2;
        } else if (dim < sizes.size()) {
            ++dim;
        } else {
            throw std::runtime_error("what the hell");
        }
    }
}

template <class Support, std::size_t N, class Functor, class... DDoms1d>
void host_for_each_block_impl(
        Support const& domain,
        std::array<ddc::DiscreteVectorElement, N> const& nb_blocks_per_dim,
        Functor const& f,
        DDoms1d const&... ddoms) noexcept
{
    static constexpr std::size_t I = sizeof...(DDoms1d);
    if constexpr (I == N) {
        f(Support(ddoms...));
    } else {
        using DDim = ddc::type_seq_element_t<I, ddc::to_type_seq_t<Support>>;
        std::size_t const block = domain.template extent<DDim>() / nb_blocks_per_dim[I];
        std::size_t const rem = domain.template extent<DDim>() - nb_blocks_per_dim[I] * block;
        ddc::DiscreteElement<DDim> front(domain.front());
        for (std::size_t ib = 0; ib < nb_blocks_per_dim[I]; ++ib) {
            ddc::DiscreteVector<DDim> const size(block + (ib < rem ? 1 : 0));
            host_for_each_block_impl(
                    domain,
                    nb_blocks_per_dim,
                    f,
                    ddoms...,
                    ddc::DiscreteDomain<DDim>(front, size));
            front += size;
        }
    }
}

template <class Support, std::size_t N, class Functor, class... DDoms1d>
void host_for_each_block(
        Support const& domain,
        std::array<ddc::DiscreteVectorElement, N> const& nb_blocks_per_dim,
        Functor const& f) noexcept
{
    host_for_each_block_impl(domain, nb_blocks_per_dim, f);
}

} // namespace detail

template <class Support, class Functor, class... DDoms1d>
void host_for_each_block(Support const& domain, std::size_t nb_blocks, Functor const& f) noexcept
{
    std::array<ddc::DiscreteVectorElement, Support::rank()> nb_blocks_per_dim {};
    detail::distribute_blocks(nb_blocks, ddc::detail::array(domain.extents()), nb_blocks_per_dim);
    detail::host_for_each_block_impl(domain, nb_blocks_per_dim, f);
}

} // namespace ddcHelper
