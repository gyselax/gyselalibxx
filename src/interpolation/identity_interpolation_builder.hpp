// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_aliases.hpp"

/**
 * @brief A builder class for copying data.
 *
 * A class which contains an operator () which can be used to build an interpolation
 * of a function. This class handles the case where no calculations are necessary and
 * the data simply needs to be copied.
 *
 * @tparam ExecSpace The Kokkos execution space on which the spline approximation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (interpolation function and splines coefficients) is stored.
 * @tparam InterpolationDDim The discrete dimension on which interpolation points are defined.
 * @tparam Basis The basis on which the interpolation is constructed.
 */
template <class ExecSpace, class MemorySpace, class InterpolationDDim, class Basis>
class IdentityInterpolationBuilder
{
public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The type of the interpolation continuous dimension (continuous dimension of interest) used by this class.
    using continuous_dimension_type = typename InterpolationDDim::continuous_dimension_type;

    /// @brief The type of the interpolation discrete dimension (discrete dimension of interest) used by this class.
    using interpolation_discrete_dimension_type = InterpolationDDim;
    /// @brief The type of the domain for the 1D interpolation mesh used by this class.
    using interpolation_domain_type = IdxRange<interpolation_discrete_dimension_type>;

    /// @brief The grid on which the interpolation coefficients should be provided.
    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    /**
     * @brief The type of the whole domain representing interpolation points.
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     */
    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_interpolation_domain_type = BatchedInterpolationGrid;

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
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batch_domain_type = ddc::
            remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;

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
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_basis_domain_type = ddc::replace_dim_of_t<
            BatchedInterpolationGrid,
            interpolation_discrete_dimension_type,
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
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_derivs_domain_type = ddc::
            remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;

public:
    /// @brief The number of equations defining the boundary condition at the lower bound.
    static constexpr int s_nbc_xmin = 0;

    /// @brief The number of equations defining the boundary condition at the upper bound.
    static constexpr int s_nbc_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

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
    template <class DataType, class Layout, class BatchedInterpolationGrid>
    void operator()(
            Field<DataType,
                  batched_basis_domain_type<BatchedInterpolationGrid>,
                  memory_space,
                  Layout> coeffs,
            ConstField<DataType, BatchedInterpolationGrid, memory_space, Layout> vals,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_domain_type<BatchedInterpolationGrid>,
                    memory_space,
                    Layout>> derivs_xmax
            = std::nullopt) const
    {
        IdxRange<basis_domain_type> bp_idx_range
                = ddc::discrete_space<Basis>().break_point_domain().remove_last(
                        IdxStep<basis_domain_type>(static_cast<int>(Basis::is_periodic())));
        Kokkos::deep_copy(
                coeffs[bp_idx_range].allocation_kokkos_view(),
                vals.allocation_kokkos_view());
        if constexpr (Basis::is_periodic()) {
            IdxRange<basis_domain_type> extended_domain(
                    ddc::discrete_space<Basis>().full_domain().remove_first(
                            bp_idx_range.extents()));
            typename BatchedInterpolationGrid::discrete_vector_type nrepeat(extended_domain.size());
            BatchedInterpolationGrid repeat_domain(get_idx_range(vals).take_first(nrepeat));
            Kokkos::deep_copy(
                    coeffs[extended_domain].allocation_kokkos_view(),
                    vals[repeat_domain].allocation_kokkos_view());
        }
    }
};
