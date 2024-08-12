// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"

/**
 * @brief A helper class which provides an operator which gives unique values for multi-D indexes
 */
template <class... Grid1D>
class GridBuilder
{
private:
    using type_seq = ddc::detail::TypeSeq<Grid1D...>;

    std::array<std::size_t, sizeof...(Grid1D)> steps;

    IdxRange<Grid1D...> m_idx_range;

public:
    /**
     * @brief A constructor for GridBuilder.
     *
     * @param[in] idx_range The valid index ranges for the operator.
     */
    GridBuilder(IdxRange<Grid1D...> idx_range) : m_idx_range(idx_range)
    {
        fill_steps<0>(idx_range);
    }

    /**
     * @brief An operator returning a unique value.
     * @param[in] idx An index identifying the requested value.
     * @return A unique value to appear in a field at the provided index.
     */
    KOKKOS_FUNCTION std::size_t operator()(Idx<Grid1D...> idx) const
    {
        return (((ddc::select<Grid1D>(idx) - ddc::select<Grid1D>(m_idx_range.front())).value()
                 * steps[ddc::type_seq_rank_v<Grid1D, type_seq>])
                + ...);
    }

private:
    template <int I, class HeadGrid1D, class... TailGrid1D>
    void fill_steps(IdxRange<HeadGrid1D, TailGrid1D...> idx_range)
    {
        if constexpr (I == (sizeof...(Grid1D) - 1)) {
            steps[I] = 1;
        } else if constexpr (I < sizeof...(Grid1D)) {
            IdxRange<TailGrid1D...> sub_idx_range(idx_range);
            steps[I] = sub_idx_range.size();
            fill_steps<I + 1>(sub_idx_range);
        }
    }
};
