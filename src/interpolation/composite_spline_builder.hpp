// SPDX-License-Identifier: MIT
#pragma once

#include <type_traits>
#include <utility>

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

/**
 * @brief A 2D spline builder which combines two 1D builders.
 *
 * The first 1D builder is used to compute the spline coefficients along the first dimension
 * (e.g. using an operator which distributes the work with MPI). The second builder (a
 * ddc::SplineBuilder) is then used to compute the spline coefficients along the second
 * dimension. The algorithm is the same as the one used in ddc::SplineBuilder2D.
 *
 * The dimension handled by the second builder must be fully available on each MPI rank.
 *
 * @tparam Builder1 The 1D builder used along the first dimension.
 * @tparam Builder2 The 1D builder (ddc::SplineBuilder) used along the second dimension.
 */
template <concepts::InterpolationBuilder1D Builder1, concepts::InterpolationBuilder1D Builder2>
class CompositeSplineBuilder2D
{
    using Traits1 = InterpolationBuilderTraits<Builder1>;
    using Traits2 = InterpolationBuilderTraits<Builder2>;

    static_assert(std::is_same_v<typename Builder1::exec_space, typename Builder2::exec_space>);
    static_assert(std::is_same_v<typename Builder1::memory_space, typename Builder2::memory_space>);
    static_assert(std::is_same_v<typename Traits1::data_type, typename Traits2::data_type>);

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = typename Builder1::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = typename Builder1::memory_space;

    /// @brief The data type that the data is saved on.
    using data_type = typename Traits1::data_type;

    /// @brief The type of the builder used along the first dimension.
    using builder_type1 = Builder1;

    /// @brief The type of the builder used along the second dimension.
    using builder_type2 = Builder2;

    /// @brief The 2D index range for the interpolation mesh.
    using interpolation_idx_range_type
            = IdxRange<interpolation_grid_t<Builder1>, interpolation_grid_t<Builder2>>;

    /// @brief The 2D index range for the interpolation coefficients.
    using coeff_idx_range_type = ddc::combine_t<
            typename Traits1::coeff_idx_range_type,
            typename Traits2::coeff_idx_range_type>;

    /**
     * @brief The batched index range on which the coefficients are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type = typename Traits2::template batched_basis_idx_range_type<
            typename Traits1::template batched_basis_idx_range_type<IdxRangeBatchedInterpolation>>;

    /**
     * @brief The batched index range on which the derivatives in the first dimension are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type1 =
            typename Traits1::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;

    /**
     * @brief The batched index range on which the derivatives in the second dimension are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type2 =
            typename Traits2::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;

    /**
     * @brief The batched index range on which the cross-derivatives are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type = typename Traits1::template batched_derivs_idx_range_type<
            typename Traits2::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>>;

private:
    Builder1 const& m_builder1;
    Builder2 const& m_builder2;

public:
    /**
     * @brief Build a CompositeSplineBuilder2D from two 1D builders.
     *
     * @param[in] builder1 The builder used along the first dimension.
     * @param[in] builder2 The builder used along the second dimension.
     */
    CompositeSplineBuilder2D(Builder1 const& builder1, Builder2 const& builder2)
        : m_builder1(builder1)
        , m_builder2(builder2)
    {
    }

    /**
     * @brief Get the batched index range on which the coefficients are defined.
     *
     * @param[in] batched_interpolation_idx_range The batched index range on which the interpolation points are defined.
     *
     * @return The batched index range on which the coefficients are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    batched_basis_idx_range_type<IdxRangeBatchedInterpolation> batched_basis_idx_range(
            IdxRangeBatchedInterpolation const& batched_interpolation_idx_range) const noexcept
    {
        return ::batched_basis_idx_range(
                m_builder2,
                ::batched_basis_idx_range(m_builder1, batched_interpolation_idx_range));
    }
};
