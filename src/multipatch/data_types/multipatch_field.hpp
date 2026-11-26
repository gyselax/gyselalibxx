// SPDX-License-Identifier: MIT

#pragma once
#include "multipatch_type.hpp"
#include "patch.hpp"
#include "vector_field.hpp"


template <class T>
inline constexpr bool enable_multipatch_field = false;

template <class T>
inline constexpr bool is_multipatch_field_v
        = enable_multipatch_field<std::remove_const_t<std::remove_reference_t<T>>>;

namespace detail {

template <class TypeSeqPatches, class TypeSeqInternalTypes>
class MultipatchField;

/**
 * @brief A class to store field objects on patches.
 *
 * On a multipatch domain when we have objects and types defined on different patches, e.g. fields.
 * They can be stored in this class and then be accessed by the patch they are defined
 * on.
 *
 * This class should not be used directly to instantiate objects. Instead the global alias
 * MultipatchFieldOf should be used.
 * I.e. use `MultipatchFieldOf<TypeOnPatch, Patches...>` instead of
 * `detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>`.
 *
 * @tparam Patches The patches on which the objects are defined.
 * @tparam InternalTypes The types of the field objects that are stored on the given patches.
 *                 The order must match that of the patches.
 *
 * @warning The objects have to be defined on different patches. Otherwise retrieving
 *          them by their patch is ill-defined.
 */
template <class... Patches, class... InternalTypes>
class MultipatchField<ddc::TypeSeq<Patches...>, ddc::TypeSeq<InternalTypes...>>
    : public MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<InternalTypes...>>
{
    static_assert(
            (has_data_access_methods_v<InternalTypes> && ...),
            "The MultipatchField type should only contain instances of objects that can be "
            "manipulated like fields.");

public:
    /// @brief The MultipatchType from which this class inherits
    using base_type = MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<InternalTypes...>>;

    /// @brief A tag storing the order of Patches in this MultipatchField
    using typename base_type::PatchOrdering;

    /// @brief The type of the object stored on the given patch.
    template <class Patch>
    using TypeOnPatch = typename base_type::template TypeOnPatch<Patch>;

    static_assert(
            !is_mem_type_v<typename base_type::example_element>,
            "For correct GPU handling a FieldMem object must be saved in a MultipatchFieldMem "
            "type.");

public:
    /// The type of a modifiable reference to this multipatch field
    using span_type
            = MultipatchField<PatchOrdering, ddc::TypeSeq<typename InternalTypes::span_type...>>;
    /// The type of a constant reference to this multipatch field
    using view_type
            = MultipatchField<PatchOrdering, ddc::TypeSeq<typename InternalTypes::view_type...>>;
    /// The type of the index ranges that can be used to access this field.
    using discrete_domain_type = MultipatchType<
            PatchOrdering,
            ddc::TypeSeq<typename InternalTypes::discrete_domain_type...>>;
    /// The memory space (CPU/GPU) where the data is saved.
    using memory_space = typename base_type::example_element::memory_space;
    /// The type of the elements inside the field.
    using element_type = typename base_type::example_element::element_type;

public:
    /**
     * Instantiate the MultipatchField class from an arbitrary number of objects.
     *
     * @param args The objects to be stored in the class.
     */
    explicit KOKKOS_FUNCTION MultipatchField(InternalTypes... args) : base_type(args...) {}

    /**
     * Create a MultipatchField class by copying an instance of another compatible MultipatchField.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects which can be converted to the correct type.
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <class MultipatchObj, std::enable_if_t<!is_mem_type_v<MultipatchObj>, bool> = true>
    explicit KOKKOS_FUNCTION MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /**
     * Create a MultipatchField class from a compatible MultipatchFieldMem.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects of the correct type.
     *
     * @param other The MultipatchFieldMem being accessed.
     */
    template <class MultipatchObj, std::enable_if_t<is_mem_type_v<MultipatchObj>, bool> = true>
    explicit MultipatchField(MultipatchObj& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(is_multipatch_type_v<MultipatchObj>);
    }

    /**
     * Create a MultipatchField class by copying an instance of another compatible MultipatchField.
     *
     * A compatible MultipatchField is one which uses all the patches used by this class. The object
     * being copied may include more patches than this MultipatchField. Further the original
     * MultipatchField must store objects of the correct type (the type template may be different
     * but return the same type depending on how it is designed.
     *
     * This function is not explicit as it is helpful to be able to change between equivalent multipatch
     * definitions if the internal type is the same but the definition comes from different locations in
     * the code.
     *
     * @param other The equivalent MultipatchField being copied.
     */
    template <class OPatchSeq, class OTypeSeq>
    KOKKOS_FUNCTION MultipatchField(MultipatchField<OPatchSeq, OTypeSeq> const& other)
        : base_type(InternalTypes(other.template get<Patches>())...)
    {
        static_assert(
                ddc::type_seq_contains_v<PatchOrdering, OPatchSeq>,
                "The type being copied does not contain all the required patches");
        static_assert(
                ddc::type_seq_contains_v<typename base_type::InternalTypeOrdering, OTypeSeq>,
                "MultipatchFields are not equivalent");
    }


    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField const&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION MultipatchField(MultipatchField&&) noexcept = default;

    KOKKOS_DEFAULTED_FUNCTION ~MultipatchField() noexcept = default;

    /**
     * Retrieve an object from the patch that it is defined on.
     *
     * @tparam Patch The patch of the object to be returned.
     * @return The object on the given patch.
     */
    template <class PatchType>
    KOKKOS_FUNCTION auto get() const requires(is_patch_v<PatchType>)
    {
        return std::get<ddc::type_seq_rank_v<PatchType, PatchOrdering>>(base_type::m_tuple);
    }

    /**
     * @brief Get a MultipatchType containing the index ranges on which the fields are defined.
     *
     * @returns The set of index ranges on which the set of fields stored in this class are defined.
     */
    KOKKOS_FUNCTION discrete_domain_type idx_range() const
    {
        return discrete_domain_type(get_idx_range(get<Patches>())...);
    }

    /**
     * @brief Get a MultipatchField containing modifiable fields.
     *
     * @returns A set of modifiable fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION span_type get_field()
    {
        return span_type(::get_field(get<Patches>())...);
    }

    /**
     * @brief Get a MultipatchField containing constant fields so the values cannot be modified.
     *
     * @returns A set of constant fields providing access to the fields stored in this class.
     */
    KOKKOS_FUNCTION view_type get_const_field() const
    {
        return view_type(::get_const_field(get<Patches>())...);
    }

    /**
     * @brief Get the Field describing the component in the QueryTag direction.
     *
     * @return The field in the specified direction.
     */
    template <class QueryTag>
    inline constexpr auto get() const noexcept
            requires((is_vector_field_v<InternalTypes> && ...) && !is_patch_v<QueryTag>)
    {
        static_assert((ddc::in_tags_v<QueryTag, typename InternalTypes::NDTypeTag> && ...));
        using FieldType = MultipatchField<
                PatchOrdering,
                ddc::TypeSeq<typename InternalTypes::chunk_span_type...>>;
        return FieldType(ddcHelper::get<QueryTag>(std::get<InternalTypes>(base_type::m_tuple))...);
    }
};

} // namespace detail

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_type<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_data_access_methods<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

template <class TypeSeqPatches, class TypeSeqInternalTypes>
inline constexpr bool enable_multipatch_field<
        detail::MultipatchField<TypeSeqPatches, TypeSeqInternalTypes>> = true;

/**
 * @brief A class to store field objects on patches.
 *
 * See detail::MultipatchField for more details. This alias creates the type T<Patch> for each of
 * the patches. Two type templates which return the same type for each patch therefore lead to the
 * same MultipatchField.
 *
 * @tparam T The type of the fields/derivative fields/vector fields that are stored on the given patches.
 * @tparam Patches The patches of the objects in the same order of the patches
 *                 that the given objects are defined on.
 */
template <template <typename P> typename T, concepts::Patch... Patches>
using MultipatchFieldOf
        = detail::MultipatchField<ddc::TypeSeq<Patches...>, ddc::TypeSeq<T<Patches>...>>;

/**
 * @brief A class to store field objects on patches.
 *
 * This alias is equivalent to MultipatchFieldOf but takes the patches as a ddc::TypeSeq.
 * I.e. `MultipatchFieldOfFromSeq<T, ddc::TypeSeq<Patches...>>` is the same type as
 * `MultipatchFieldOf<T, Patches...>`.
 *
 * @tparam T The type of the fields/derivative fields/vector fields that are stored on the given patches.
 * @tparam PatchSeq A ddc::TypeSeq of the patches of the objects.
 */
template <template <typename P> typename T, ddc::concepts::type_seq PatchSeq>
using MultipatchFieldOfFromSeq =
        typename detail::MultipatchFromPatchSeq<detail::MultipatchField, T, PatchSeq>::type;

namespace detail {

/**
 * @brief A helper to build a MultipatchField of Fields from a MultipatchType of index ranges.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <class ElementType, class MultipatchIdxRange, class MemorySpace, class LayoutStridedPolicy>
struct MultipatchFieldFromIdxRange;

/// @brief Specialisation of MultipatchFieldFromIdxRange to unpack the patches and index ranges.
template <
        class ElementType,
        class... Patches,
        class... IdxRanges,
        class MemorySpace,
        class LayoutStridedPolicy>
struct MultipatchFieldFromIdxRange<
        ElementType,
        MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<IdxRanges...>>,
        MemorySpace,
        LayoutStridedPolicy>
{
    static_assert(
            (ddc::is_discrete_domain_v<IdxRanges> && ...),
            "The MultipatchType must contain an index range on each patch.");
    /// The MultipatchField of Fields defined on the index ranges.
    using type = MultipatchField<
            ddc::TypeSeq<Patches...>,
            ddc::TypeSeq<Field<ElementType, IdxRanges, MemorySpace, LayoutStridedPolicy>...>>;
};

} // namespace detail

/**
 * @brief A class to store a Field on each patch.
 *
 * This alias is templated in the same way as Field, with the index range replaced by a
 * MultipatchType of index ranges (e.g. `MultipatchType<IdxRangeOnPatch, Patches...>`). The
 * resulting type is the same as `MultipatchFieldOf<T, Patches...>` where
 * `T<Patch> = Field<ElementType, IdxRangeOnPatch<Patch>, MemorySpace, LayoutStridedPolicy>`.
 *
 * VectorFields should be stored using MultipatchVectorField. Other fields which are not plain
 * Fields (e.g. DerivField) should be stored using MultipatchFieldOf.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class ElementType,
        class MultipatchIdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using MultipatchField = typename detail::MultipatchFieldFromIdxRange<
        ElementType,
        MultipatchIdxRange,
        MemorySpace,
        LayoutStridedPolicy>::type;

/**
 * @brief A class to store a Field of doubles on each patch.
 * See MultipatchField.
 *
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class MultipatchIdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DMultipatchField
        = MultipatchField<double, MultipatchIdxRange, MemorySpace, LayoutStridedPolicy>;

/**
 * @brief A class to store a constant Field on each patch.
 * See MultipatchField.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class ElementType,
        class MultipatchIdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using MultipatchConstField
        = MultipatchField<ElementType const, MultipatchIdxRange, MemorySpace, LayoutStridedPolicy>;

/**
 * @brief A class to store a constant Field of doubles on each patch.
 * See MultipatchField.
 *
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class MultipatchIdxRange,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DMultipatchConstField
        = MultipatchConstField<double, MultipatchIdxRange, MemorySpace, LayoutStridedPolicy>;

namespace detail {

/**
 * @brief A helper to build a MultipatchField of VectorFields from a MultipatchType of index ranges.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam VectorIndexSetType A VectorIndexSet describing the components of the vector field.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class ElementType,
        class MultipatchIdxRange,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
struct MultipatchVectorFieldFromIdxRange;

/// @brief Specialisation of MultipatchVectorFieldFromIdxRange to unpack the patches and index ranges.
template <
        class ElementType,
        class... Patches,
        class... IdxRanges,
        class VectorIndexSetType,
        class MemorySpace,
        class LayoutStridedPolicy>
struct MultipatchVectorFieldFromIdxRange<
        ElementType,
        MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<IdxRanges...>>,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>
{
    static_assert(
            (ddc::is_discrete_domain_v<IdxRanges> && ...),
            "The MultipatchType must contain an index range on each patch.");
    /// The MultipatchField of VectorFields defined on the index ranges.
    using type = MultipatchField<
            ddc::TypeSeq<Patches...>,
            ddc::TypeSeq<VectorField<
                    ElementType,
                    IdxRanges,
                    VectorIndexSetType,
                    MemorySpace,
                    LayoutStridedPolicy>...>>;
};

} // namespace detail

/**
 * @brief A class to store a VectorField on each patch.
 *
 * This alias is templated in the same way as VectorField, with the index range replaced by a
 * MultipatchType of index ranges (e.g. `MultipatchType<IdxRangeOnPatch, Patches...>`). The
 * resulting type is the same as `MultipatchFieldOf<T, Patches...>` where
 * `T<Patch> = VectorField<ElementType, IdxRangeOnPatch<Patch>, VectorIndexSetType, MemorySpace,
 * LayoutStridedPolicy>`.
 *
 * The same components (VectorIndexSetType) are used on every patch. Vector fields whose
 * components differ from one patch to another should be stored using MultipatchFieldOf.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam VectorIndexSetType A VectorIndexSet describing the components of the vector field.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class ElementType,
        class MultipatchIdxRange,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using MultipatchVectorField = typename detail::MultipatchVectorFieldFromIdxRange<
        ElementType,
        MultipatchIdxRange,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>::type;

/**
 * @brief A class to store a VectorField of doubles on each patch.
 * See MultipatchVectorField.
 *
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam VectorIndexSetType A VectorIndexSet describing the components of the vector field.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class MultipatchIdxRange,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DMultipatchVectorField = MultipatchVectorField<
        double,
        MultipatchIdxRange,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>;

/**
 * @brief A class to store a constant VectorField on each patch.
 * See MultipatchVectorField.
 *
 * @tparam ElementType The type of the elements inside the fields.
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam VectorIndexSetType A VectorIndexSet describing the components of the vector field.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class ElementType,
        class MultipatchIdxRange,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using MultipatchVectorConstField = MultipatchVectorField<
        ElementType const,
        MultipatchIdxRange,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>;

/**
 * @brief A class to store a constant VectorField of doubles on each patch.
 * See MultipatchVectorField.
 *
 * @tparam MultipatchIdxRange A MultipatchType containing the index range on each patch.
 * @tparam VectorIndexSetType A VectorIndexSet describing the components of the vector field.
 * @tparam MemorySpace The memory space (CPU/GPU) where the data is saved.
 * @tparam LayoutStridedPolicy The layout of the data in memory.
 */
template <
        class MultipatchIdxRange,
        class VectorIndexSetType,
        class MemorySpace = Kokkos::DefaultExecutionSpace::memory_space,
        class LayoutStridedPolicy = Kokkos::layout_right>
using DMultipatchVectorConstField = MultipatchVectorConstField<
        double,
        MultipatchIdxRange,
        VectorIndexSetType,
        MemorySpace,
        LayoutStridedPolicy>;

template <class... VectorFieldType, class PatchTypeSeq>
inline constexpr bool enable_vector_field<detail::MultipatchField<
        PatchTypeSeq,
        ddc::TypeSeq<VectorFieldType...>>> = (is_vector_field_v<VectorFieldType> && ...);

namespace ddcHelper {

/**
 * @brief Copy the data from one MultipatchField into another
 * @param dst The MultipatchField that the data will be copied to.
 * @param src The MultipatchField that the data will be copied from.
 */
template <class... Patches, class TypeSeqInternalTypes1, class TypeSeqInternalTypes2>
void deepcopy(
        detail::MultipatchField<ddc::TypeSeq<Patches...>, TypeSeqInternalTypes1> dst,
        detail::MultipatchField<ddc::TypeSeq<Patches...>, TypeSeqInternalTypes2> src)
{
    using DstType = detail::MultipatchField<ddc::TypeSeq<Patches...>, TypeSeqInternalTypes1>;
    if constexpr (ddc::is_chunk_v<typename DstType::example_element>) {
        (ddc::parallel_deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    } else {
        (ddcHelper::deepcopy(dst.template get<Patches>(), src.template get<Patches>()), ...);
    }
}

} // namespace ddcHelper
