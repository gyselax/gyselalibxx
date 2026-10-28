// SPDX-License-Identifier: MIT
#pragma once

#include <optional>
#include <type_traits>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"


/**
 * @brief A 3D spline builder which combines a 1D builder and a 2D builder.
 *
 * The 1D builder is used to compute the spline coefficients along the first dimension
 * (e.g. using an operator which distributes the work with MPI). The 2D builder (a
 * ddc::SplineBuilder2D) is then used to compute the spline coefficients along the second
 * and third dimensions. The algorithm is the same as the one used in ddc::SplineBuilder3D.
 *
 * The dimensions handled by the 2D builder must be fully available on each MPI rank.
 *
 * @tparam Builder1 The 1D builder used along the first dimension.
 * @tparam Builder23 The 2D builder (ddc::SplineBuilder2D) used along the second and third dimensions.
 */
template <concepts::InterpolationBuilder1D Builder1, concepts::InterpolationBuilder Builder23>
class CompositeSplineBuilder3D
{
    using Traits1 = InterpolationBuilderTraits<Builder1>;
    using Traits23 = InterpolationBuilderTraits<Builder23>;
    // The per-dimension properties of the 2D builder are obtained from the 1D builders that it uses
    using Builder2 = typename Builder23::builder_type1;
    using Builder3 = typename Builder23::builder_type2;
    using Traits2 = InterpolationBuilderTraits<Builder2>;
    using Traits3 = InterpolationBuilderTraits<Builder3>;

    static_assert(Traits23::rank() == 2);
    static_assert(std::is_same_v<typename Builder1::exec_space, typename Builder23::exec_space>);
    static_assert(std::is_same_v<typename Builder1::memory_space, typename Builder23::memory_space>);
    static_assert(std::is_same_v<typename Traits1::data_type, typename Traits23::data_type>);

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = typename Builder1::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = typename Builder1::memory_space;

    /// @brief The data type that the data is saved on.
    using data_type = typename Traits1::data_type;

    /// @brief The type of the builder used along the first dimension.
    using builder_type1 = Builder1;

    /// @brief The type of the builder used along the second and third dimensions.
    using builder_type_2_3 = Builder23;

    /// @brief The 3D index range for the interpolation mesh.
    using interpolation_idx_range_type = ddc::combine_t<
            typename Traits1::interpolation_idx_range_type,
            typename Traits23::interpolation_idx_range_type>;

    /// @brief The 3D index range for the interpolation coefficients.
    using coeff_idx_range_type = ddc::combine_t<
            typename Traits1::coeff_idx_range_type,
            typename Traits23::coeff_idx_range_type>;

    /**
     * @brief The batched index range on which the coefficients are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type = typename Traits23::template batched_basis_idx_range_type<
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
     * @brief The batched index range on which the derivatives in the third dimension are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type3 =
            typename Traits3::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;

    /**
     * @brief The batched index range on which the cross-derivatives in the first and second
     * dimensions are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type1_2 = typename Traits1::template batched_derivs_idx_range_type<
            batched_derivs_idx_range_type2<IdxRangeBatchedInterpolation>>;

    /**
     * @brief The batched index range on which the cross-derivatives in the second and third
     * dimensions are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type2_3 = typename Traits2::template batched_derivs_idx_range_type<
            batched_derivs_idx_range_type3<IdxRangeBatchedInterpolation>>;

    /**
     * @brief The batched index range on which the cross-derivatives in the first and third
     * dimensions are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type1_3 = typename Traits1::template batched_derivs_idx_range_type<
            batched_derivs_idx_range_type3<IdxRangeBatchedInterpolation>>;

    /**
     * @brief The batched index range on which the cross-derivatives in all three dimensions
     * are defined.
     *
     * @tparam IdxRangeBatchedInterpolation The batched index range on which the interpolation points are defined.
     */
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type = typename Traits1::template batched_derivs_idx_range_type<
            batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>>;

private:
    Builder1 const& m_builder1;
    Builder23 const& m_builder23;

public:
    /**
     * @brief Build a CompositeSplineBuilder3D from a 1D builder and a 2D builder.
     *
     * @param[in] builder1 The builder used along the first dimension.
     * @param[in] builder23 The builder used along the second and third dimensions.
     */
    CompositeSplineBuilder3D(Builder1 const& builder1, Builder23 const& builder23)
        : m_builder1(builder1)
        , m_builder23(builder23)
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
                m_builder23,
                ::batched_basis_idx_range(m_builder1, batched_interpolation_idx_range));
    }

    /**
     * @brief Compute a 3D spline approximation of a function.
     *
     * Use the values of a function and the derivatives of the function at the boundaries
     * (in the case of Hermite boundary conditions only) to calculate the coefficients of a
     * 3D spline approximation of this function.
     *
     * @param[out] coeffs The coefficients of the spline.
     * @param[in] vals The values of the function at the interpolation mesh.
     * @param[in] derivs_min1 The values of the derivatives at the lower boundary in the first dimension.
     * @param[in] derivs_max1 The values of the derivatives at the upper boundary in the first dimension.
     * @param[in] derivs_min2 The values of the derivatives at the lower boundary in the second dimension.
     * @param[in] derivs_max2 The values of the derivatives at the upper boundary in the second dimension.
     * @param[in] derivs_min3 The values of the derivatives at the lower boundary in the third dimension.
     * @param[in] derivs_max3 The values of the derivatives at the upper boundary in the third dimension.
     * @param[in] mixed_derivs_min1_min2 The values of the cross-derivatives at the lower boundary
     *      in the first dimension and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_max1_min2 The values of the cross-derivatives at the upper boundary
     *      in the first dimension and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_min1_max2 The values of the cross-derivatives at the lower boundary
     *      in the first dimension and the upper boundary in the second dimension.
     * @param[in] mixed_derivs_max1_max2 The values of the cross-derivatives at the upper boundary
     *      in the first dimension and the upper boundary in the second dimension.
     * @param[in] mixed_derivs_min2_min3 The values of the cross-derivatives at the lower boundary
     *      in the second dimension and the lower boundary in the third dimension.
     * @param[in] mixed_derivs_max2_min3 The values of the cross-derivatives at the upper boundary
     *      in the second dimension and the lower boundary in the third dimension.
     * @param[in] mixed_derivs_min2_max3 The values of the cross-derivatives at the lower boundary
     *      in the second dimension and the upper boundary in the third dimension.
     * @param[in] mixed_derivs_max2_max3 The values of the cross-derivatives at the upper boundary
     *      in the second dimension and the upper boundary in the third dimension.
     * @param[in] mixed_derivs_min1_min3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension and the lower boundary in the third dimension.
     * @param[in] mixed_derivs_max1_min3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension and the lower boundary in the third dimension.
     * @param[in] mixed_derivs_min1_max3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension and the upper boundary in the third dimension.
     * @param[in] mixed_derivs_max1_max3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension and the upper boundary in the third dimension.
     * @param[in] mixed_derivs_min1_min2_min3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension, the lower boundary in the second dimension and the lower
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_max1_min2_min3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension, the lower boundary in the second dimension and the lower
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_min1_max2_min3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension, the upper boundary in the second dimension and the lower
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_max1_max2_min3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension, the upper boundary in the second dimension and the lower
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_min1_min2_max3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension, the lower boundary in the second dimension and the upper
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_max1_min2_max3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension, the lower boundary in the second dimension and the upper
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_min1_max2_max3 The values of the cross-derivatives at the lower boundary
     *      in the first dimension, the upper boundary in the second dimension and the upper
     *      boundary in the third dimension.
     * @param[in] mixed_derivs_max1_max2_max3 The values of the cross-derivatives at the upper boundary
     *      in the first dimension, the upper boundary in the second dimension and the upper
     *      boundary in the third dimension.
     */
    template <class IdxRangeBatchedInterpolation, class Layout>
    void operator()(
            Field<data_type,
                  batched_basis_idx_range_type<IdxRangeBatchedInterpolation>,
                  memory_space,
                  Layout> coeffs,
            ConstField<data_type, IdxRangeBatchedInterpolation, memory_space, Layout> vals,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_min1
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_max1
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_min2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_max2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> derivs_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_min2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_min2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_max2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_2<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_max2
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min2_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max2_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type1_3<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_min2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_min2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_max2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_max2_min3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_min2_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_min2_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_min1_max2_max3
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    Layout>> mixed_derivs_max1_max2_max3
            = std::nullopt) const
    {
        using IdxRangeCoeffs1Derivs2 = typename Traits1::template batched_basis_idx_range_type<
                batched_derivs_idx_range_type2<IdxRangeBatchedInterpolation>>;
        using IdxRangeCoeffs1Derivs3 = typename Traits1::template batched_basis_idx_range_type<
                batched_derivs_idx_range_type3<IdxRangeBatchedInterpolation>>;
        using IdxRangeCoeffs1Derivs2_3 = typename Traits1::template batched_basis_idx_range_type<
                batched_derivs_idx_range_type2_3<IdxRangeBatchedInterpolation>>;

        constexpr bool hermite_min2 = Builder2::s_sbc_xmin == ddc::SplineBuilderClosure::HERMITE;
        constexpr bool hermite_max2 = Builder2::s_sbc_xmax == ddc::SplineBuilderClosure::HERMITE;
        constexpr bool hermite_min3 = Builder3::s_sbc_xmin == ddc::SplineBuilderClosure::HERMITE;
        constexpr bool hermite_max3 = Builder3::s_sbc_xmax == ddc::SplineBuilderClosure::HERMITE;

        // Do not return early if the batch index range is empty as Builder1 may require
        // the participation of all MPI ranks.

        // Build the coefficients along the first dimension of the derivatives at the lower
        // boundary of the second dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2, memory_space>>
                coeffs1_derivs_min2_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2, memory_space>>
                coeffs1_derivs_min2;
        if constexpr (hermite_min2) {
            coeffs1_derivs_min2_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*derivs_min2)));
            m_builder1(
                    get_field(*coeffs1_derivs_min2_alloc),
                    *derivs_min2,
                    mixed_derivs_min1_min2,
                    mixed_derivs_max1_min2);
            coeffs1_derivs_min2 = get_const_field(*coeffs1_derivs_min2_alloc);
        }

        // Build the coefficients along the first dimension of the derivatives at the upper
        // boundary of the second dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2, memory_space>>
                coeffs1_derivs_max2_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2, memory_space>>
                coeffs1_derivs_max2;
        if constexpr (hermite_max2) {
            coeffs1_derivs_max2_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*derivs_max2)));
            m_builder1(
                    get_field(*coeffs1_derivs_max2_alloc),
                    *derivs_max2,
                    mixed_derivs_min1_max2,
                    mixed_derivs_max1_max2);
            coeffs1_derivs_max2 = get_const_field(*coeffs1_derivs_max2_alloc);
        }

        // Build the coefficients along the first dimension of the derivatives at the lower
        // boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs3, memory_space>>
                coeffs1_derivs_min3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs3, memory_space>>
                coeffs1_derivs_min3;
        if constexpr (hermite_min3) {
            coeffs1_derivs_min3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*derivs_min3)));
            m_builder1(
                    get_field(*coeffs1_derivs_min3_alloc),
                    *derivs_min3,
                    mixed_derivs_min1_min3,
                    mixed_derivs_max1_min3);
            coeffs1_derivs_min3 = get_const_field(*coeffs1_derivs_min3_alloc);
        }

        // Build the coefficients along the first dimension of the derivatives at the upper
        // boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs3, memory_space>>
                coeffs1_derivs_max3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs3, memory_space>>
                coeffs1_derivs_max3;
        if constexpr (hermite_max3) {
            coeffs1_derivs_max3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*derivs_max3)));
            m_builder1(
                    get_field(*coeffs1_derivs_max3_alloc),
                    *derivs_max3,
                    mixed_derivs_min1_max3,
                    mixed_derivs_max1_max3);
            coeffs1_derivs_max3 = get_const_field(*coeffs1_derivs_max3_alloc);
        }

        // The cross-derivatives in the second and third dimensions are only used by Builder23
        // if both of these dimensions use Hermite boundary conditions. They must therefore
        // only be accessed if both dimensions are Hermite (&&, not ||). This differs from
        // ddc::SplineBuilder3D which requires them as soon as one of the dimensions is Hermite
        // and fails if only the third dimension is Hermite.

        // Build the coefficients along the first dimension of the cross-derivatives at the
        // lower boundary of the second dimension and the lower boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_min2_min3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_min2_min3;
        if constexpr (hermite_min2 && hermite_min3) {
            coeffs1_derivs_min2_min3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*mixed_derivs_min2_min3)));
            m_builder1(
                    get_field(*coeffs1_derivs_min2_min3_alloc),
                    *mixed_derivs_min2_min3,
                    mixed_derivs_min1_min2_min3,
                    mixed_derivs_max1_min2_min3);
            coeffs1_derivs_min2_min3 = get_const_field(*coeffs1_derivs_min2_min3_alloc);
        }

        // Build the coefficients along the first dimension of the cross-derivatives at the
        // upper boundary of the second dimension and the lower boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_max2_min3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_max2_min3;
        if constexpr (hermite_max2 && hermite_min3) {
            coeffs1_derivs_max2_min3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*mixed_derivs_max2_min3)));
            m_builder1(
                    get_field(*coeffs1_derivs_max2_min3_alloc),
                    *mixed_derivs_max2_min3,
                    mixed_derivs_min1_max2_min3,
                    mixed_derivs_max1_max2_min3);
            coeffs1_derivs_max2_min3 = get_const_field(*coeffs1_derivs_max2_min3_alloc);
        }

        // Build the coefficients along the first dimension of the cross-derivatives at the
        // lower boundary of the second dimension and the upper boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_min2_max3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_min2_max3;
        if constexpr (hermite_min2 && hermite_max3) {
            coeffs1_derivs_min2_max3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*mixed_derivs_min2_max3)));
            m_builder1(
                    get_field(*coeffs1_derivs_min2_max3_alloc),
                    *mixed_derivs_min2_max3,
                    mixed_derivs_min1_min2_max3,
                    mixed_derivs_max1_min2_max3);
            coeffs1_derivs_min2_max3 = get_const_field(*coeffs1_derivs_min2_max3_alloc);
        }

        // Build the coefficients along the first dimension of the cross-derivatives at the
        // upper boundary of the second dimension and the upper boundary of the third dimension
        std::optional<FieldMem<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_max2_max3_alloc;
        std::optional<ConstField<data_type, IdxRangeCoeffs1Derivs2_3, memory_space>>
                coeffs1_derivs_max2_max3;
        if constexpr (hermite_max2 && hermite_max3) {
            coeffs1_derivs_max2_max3_alloc.emplace(
                    ::batched_basis_idx_range(m_builder1, get_idx_range(*mixed_derivs_max2_max3)));
            m_builder1(
                    get_field(*coeffs1_derivs_max2_max3_alloc),
                    *mixed_derivs_max2_max3,
                    mixed_derivs_min1_max2_max3,
                    mixed_derivs_max1_max2_max3);
            coeffs1_derivs_max2_max3 = get_const_field(*coeffs1_derivs_max2_max3_alloc);
        }

        // Build the coefficients along the first dimension
        FieldMem<data_type,
                 typename Traits1::template batched_basis_idx_range_type<
                         IdxRangeBatchedInterpolation>,
                 memory_space>
                coeffs1_alloc(::batched_basis_idx_range(m_builder1, get_idx_range(vals)));
        m_builder1(get_field(coeffs1_alloc), vals, derivs_min1, derivs_max1);

        // Build the coefficients along the second and third dimensions
        m_builder23(
                coeffs,
                get_const_field(coeffs1_alloc),
                coeffs1_derivs_min2,
                coeffs1_derivs_max2,
                coeffs1_derivs_min3,
                coeffs1_derivs_max3,
                coeffs1_derivs_min2_min3,
                coeffs1_derivs_max2_min3,
                coeffs1_derivs_min2_max3,
                coeffs1_derivs_max2_max3);
    }
};
