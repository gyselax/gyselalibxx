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
    using interpolation_domain_type = IdxRange<interpolation_grid_type>;

    /// @brief The grid on which the interpolation coefficients should be provided.
    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    /**
     * @brief The type of the whole domain representing interpolation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange_,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange_>>>
    using batched_interpolation_domain_type = BatchedInterpolationIdxRange_;

    /**
     * @brief The type of the batch domain (obtained by removing the dimension of interest
     * from the whole domain).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     *
     * Example: For batched_interpolation_domain_type = DiscreteDomain<X,Y,Z> and a dimension of interest Y,
     * this is DiscreteDomain<X,Z>
     */
    template <
            class BatchedInterpolationIdxRange_,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange_>>>
    using batch_domain_type
            = ddc::remove_dims_of_t<BatchedInterpolationIdxRange_, interpolation_grid_type>;

    /**
     * @brief The type of the whole interpolation domain (cartesian product of 1D interpolation domain
     * and batch domain) preserving the underlying memory layout (order of dimensions).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     *
     * Example: For batched_interpolation_domain_type = DiscreteDomain<X,Y,Z> and a dimension of interest Y
     * (associated to a Basis tag LagrangeBasisY), this is DiscreteDomain<X,LagrangeBasisY,Z>.
     */
    template <
            class BatchedInterpolationIdxRange_,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange_>>>
    using batched_basis_domain_type = ddc::replace_dim_of_t<
            BatchedInterpolationIdxRange_,
            interpolation_grid_type,
            basis_domain_type>;

    /**
     * @brief The typeof the derivatives
     *
     * The type of the derivatives that need to be provided to this method. No derivatives
     * are required but this is included for interoperability with other interpolation
     * builder classes.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationIdxRange_,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationIdxRange_>>>
    using batched_derivs_domain_type
            = ddc::remove_dims_of_t<BatchedInterpolationIdxRange_, interpolation_grid_type>;

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
            Field<DataType, batched_basis_domain_type<BatchedInterpolationIdxRange>, memory_space>
                    coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmax
            = std::nullopt) const = 0;
};
