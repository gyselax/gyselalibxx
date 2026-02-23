// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"

template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class InterpolationGrid,
        class Basis,
        class BatchedInterpolationIdxRange>
class IInterpolationBuilder
{
public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The type of the interpolation continuous dimension (continuous dimension of interest) used by this class.
    using continuous_dimension_type = typename InterpolationGrid::continuous_dimension_type;

    /// @brief The type of the interpolation discrete dimension (discrete dimension of interest) used by this class.
    using interpolation_grid_type = InterpolationGrid;

    /// @brief The type of the domain for the 1D interpolation mesh used by this class.
    using interpolation_idx_range_type = IdxRange<interpolation_grid_type>;

    /// @brief The grid on which the interpolation coefficients should be provided.
    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    /// @brief The type of the Deriv dimension at the boundaries.
    using deriv_type = ddc::Deriv<continuous_dimension_type>;

    /**
     * @brief The type of the whole domain representing interpolation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    using batched_interpolation_idx_range_type = BatchedInterpolationIdxRange;

    /**
     * @brief The type of the batch domain (obtained by removing the dimension of interest
     * from the whole domain).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     *
     * Example: For batched_interpolation_idx_range_type = DiscreteDomain<X,Y,Z> and a dimension of interest Y,
     * this is DiscreteDomain<X,Z>
     */
    using batch_domain_type
            = ddc::remove_dims_of_t<BatchedInterpolationIdxRange, interpolation_grid_type>;

    /**
     * @brief The type of the whole interpolation domain (cartesian product of 1D interpolation domain
     * and batch domain) preserving the underlying memory layout (order of dimensions).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     *
     * Example: For batched_interpolation_idx_range_type = DiscreteDomain<X,Y,Z> and a dimension of interest Y
     * (associated to a Basis tag LagrangeBasisY), this is DiscreteDomain<X,LagrangeBasisY,Z>.
     */
    using batched_basis_idx_range_type = ddc::replace_dim_of_t<
            BatchedInterpolationIdxRange,
            interpolation_grid_type,
            basis_domain_type>;

    /**
     * @brief The type of the derivatives
     *
     * The type of the derivatives that need to be provided to this method. No derivatives
     * are required but this is included for interoperability with other interpolation
     * builder classes.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    using batched_derivs_idx_range_type = ddc::
            replace_dim_of_t<BatchedInterpolationIdxRange, interpolation_grid_type, deriv_type>;

public:
    IInterpolationBuilder() = default;

    /**
     * @brief Compute the interpolation coefficients for a function.
     *
     * No calculations are necessary and the data simply needs to be copied
     * from vals to coeffs.
     *
     * @param[out] coeffs The coefficients of the spline computed by this SplineBuilder.
     * @param[in] vals The values of the function on the interpolation mesh.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     *                  (unused in this class).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     *                  (unused in this class).
     */
    virtual void operator()(
            Field<DataType, batched_basis_idx_range_type, memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<DataType, batched_derivs_idx_range_type, memory_space>>
                    derivs_xmin
            = std::nullopt,
            std::optional<ConstField<DataType, batched_derivs_idx_range_type, memory_space>>
                    derivs_xmax
            = std::nullopt) const = 0;

    /**
     * @brief Get the whole domain on which derivatives on lower boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    virtual batched_derivs_idx_range_type batched_derivs_xmin_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept = 0;

    /**
     * @brief Get the whole domain on which derivatives on upper boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    virtual batched_derivs_idx_range_type batched_derivs_xmax_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept = 0;
};
