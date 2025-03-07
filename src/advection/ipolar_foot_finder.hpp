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
 * @tparam AdvectionDim1 The first dimension of the advection field vector.
 * @tparam AdvectionDim2 The second dimension of the advection field vector.
 * @tparam MemorySpace The memory space where the data is saved (CPU/GPU).
 */
template <
        class GridRadial,
        class GridPoloidal,
        class AdvectionDim1,
        class AdvectionDim2,
        class MemorySpace>
class IPolarFootFinder
{
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
    using X = AdvectionDim1;
    /// The continuous poloidal dimension.
    using Y = AdvectionDim2;

    /// The type of the memory space where the field is saved (CPU vs GPU).
    using memory_space = MemorySpace;

    /// The type of the index range over which the operator works.
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

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
            Field<Coord<R, Theta>, IdxRangeRTheta, memory_space> feet,
            DVectorConstField<IdxRangeRTheta, VectorIndexSet<X, Y>, memory_space> advection_field,
            double dt) const = 0;
};
