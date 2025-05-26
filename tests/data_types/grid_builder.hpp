// SPDX-License-Identifier: MIT
#pragma once
#include <cstddef>

#include <Kokkos_Core.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief A helper class which provides an operator which gives unique values for multi-D indexes
 */
template <class... Grid1D>
class GridBuilder
{
private:
    using type_seq = ddc::detail::TypeSeq<Grid1D...>;

    IdxRange<Grid1D...> m_idx_range;

    Kokkos::layout_right::mapping<Kokkos::dextents<std::size_t, sizeof...(Grid1D)>> m_mapping;

public:
    /**
     * @brief A constructor for GridBuilder.
     *
     * @param[in] idx_range The valid index ranges for the operator.
     */
    explicit GridBuilder(IdxRange<Grid1D...> idx_range)
        : m_idx_range(idx_range)
        , m_mapping(Kokkos::dextents<std::size_t, sizeof...(Grid1D)>(
                  idx_range.template extent<Grid1D>()...))
    {
    }

    /**
     * @brief An operator returning a unique value.
     * @param[in] idx An index identifying the requested value.
     * @return A unique value to appear in a field at the provided index.
     */
    KOKKOS_FUNCTION std::size_t operator()(Idx<Grid1D...> idx) const
    {
        return m_mapping(
                (ddc::select<Grid1D>(idx) - ddc::select<Grid1D>(m_idx_range.front())).value()...);
    }
};
