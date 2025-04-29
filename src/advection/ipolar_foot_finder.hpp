// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include "ddc_aliases.hpp"
#include "vector_field.hpp"

/**
 * @brief Define a base class for all the time integration methods used to find
 * the foot of a characteristic on a polar domain (a polar domain is a domain
 * defined on the @f$ (r,\theta) @f$ plane).
 *
)
 *
 * @tparam GridRadial The radial grid on which the distribution function is defined.
 * @tparam GridPoloidal The poloidial grid on which the distribution function is defined.
 * @tparam VectorIndexSetAdvDims A vector index set containing the set of dimensions
 *                              that can be used to index the advection dimensions.
 * @tparam IdxRangeBatched The index range on which this batched operator operates.
 *                              I.e. the index range of the distribution function.
 * @tparam MemorySpace The memory space where the data is saved (CPU/GPU).
 */
template <
        class GridRadial,
        class GridPoloidal,
        class VectorIndexSetAdvDims,
        class IdxRangeBatched,
        class MemorySpace>
class IPolarFootFinder
{
    // Check types
    static_assert(
            (ddc::is_non_uniform_point_sampling_v<GridRadial>)
            || (ddc::is_uniform_point_sampling_v<GridRadial>));
    static_assert(
            (ddc::is_non_uniform_point_sampling_v<GridPoloidal>)
            || (ddc::is_uniform_point_sampling_v<GridPoloidal>));
    static_assert(is_vector_index_set_v<VectorIndexSetAdvDims>);
    static_assert(ddc::is_discrete_domain_v<IdxRangeBatched>);
    static_assert(Kokkos::is_memory_space_v<MemorySpace>);

    // Check that grids make sense
    static_assert(
            ddc::in_tags_v<GridRadial, ddc::to_type_seq_t<IdxRangeBatched>>,
            "The radial grid must be found in the batched index range");
    static_assert(
            ddc::in_tags_v<GridPoloidal, ddc::to_type_seq_t<IdxRangeBatched>>,
            "The poloidal grid must be found in the batched index range");

    // Check that VectorIndexSetAdvDims makes sense
    static_assert(ddc::type_seq_size_v<VectorIndexSetAdvDims> == 2);

protected:
    /// The continuous radial dimension.
    using GridR = GridRadial;
    /// The continuous poloidal dimension.
    using GridTheta = GridPoloidal;

    /// The continuous radial dimension.
    using R = typename GridR::continuous_dimension_type;
    /// The continuous poloidal dimension.
    using Theta = typename GridTheta::continuous_dimension_type;

    /// The continuous radial dimension.
    using VectorIndexSetAdvectionDims = VectorIndexSetAdvDims;

public:
    /// The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = MemorySpace;

    /// The type of the index range over which the operator works.
    using IdxRangeOperator = IdxRangeBatched;

public:
    virtual ~IPolarFootFinder() = default;

    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * @param[in, out] feet
     *      On input: the mesh points.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the physical domain.
     * @param[in] dt
     *      The time step.
     */
    virtual void operator()(
            Field<Coord<R, Theta>, IdxRangeOperator, memory_space> feet,
            DVectorConstField<IdxRangeOperator, VectorIndexSetAdvectionDims, memory_space>
                    advection_field,
            double dt) const = 0;
};
