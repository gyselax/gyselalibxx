#pragma once

#include <array>

namespace detail {

template <class>
struct ArrayApplyer;
template <size_t... IDXS>
struct ArrayApplyer<std::index_sequence<IDXS...>>
{
    static constexpr size_t NB_ELEM = sizeof...(IDXS);

    template <typename Functor, typename Element>
    static std::array<std::invoke_result_t<Functor, const Element&>, NB_ELEM> array_apply(
            Functor&& func,
            const std::array<Element, NB_ELEM>& params)
    {
        return std::array<std::invoke_result_t<Functor, const Element&>, NB_ELEM> {
                std::forward<Functor>(func)(params[IDXS])...};
    }
};

} // namespace detail

template <typename Functor, typename Element, size_t NB_ELEM>
std::array<std::invoke_result_t<Functor, const Element&>, NB_ELEM> array_apply(
        Functor&& func,
        const std::array<Element, NB_ELEM>& params)
{
    return detail::ArrayApplyer<
            std::make_index_sequence<NB_ELEM>>::array_apply(std::forward<Functor>(func), params);
}

using RCoord = double;

template <int NDIMS>
using RCoordND = std::array<RCoord, NDIMS>;

using RCoord1D = RCoordND<1>;

using RCoord2D = RCoordND<2>;

using RCoord3D = RCoordND<3>;

class RDimension
{
    /** Periodicity of the dimension
     *
     * - in case the dimension is periodic, m_periodicity != 0.
     *   && normalize(x+n.m_periodicity)== x,
     * - otherwise, m_periodicity==0. and is unused
     */
    RCoord1D m_periodicity;

public:
    RCoord1D normalize(RCoord1D) const;
};

struct RDomain
{
    RCoord start;

    RCoord end;
};

template <int NDIMS>
using RDomainND = std::array<RDomain, NDIMS>;

using RDomain1D = RDomainND<1>;

using RDomain2D = RDomainND<2>;

using RDomain3D = RDomainND<3>;
