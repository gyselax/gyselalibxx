// SPDX-License-Identifier: MIT

#pragma once

/**
 * @brief A builder class for copying data.
 *
 * A class which contains an operator () which can be used to build a spline approximation
 * of a function. A spline approximation is represented by coefficients stored in a Chunk
 * of B-splines. The spline is constructed such that it respects the boundary conditions
 * BcLower and BcUpper, and it interpolates the function at the points on the interpolation_discrete_dimension
 * associated with interpolation_discrete_dimension_type.
 * @tparam ExecSpace The Kokkos execution space on which the spline approximation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (interpolation function and splines coefficients) is stored.
 * @tparam InterpolationDDim The discrete dimension on which interpolation points are defined.
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
    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_discrete_dimension_type>;

    using basis_domain_type = Basis::Impl<Basis, MemorySpace>::knot_grid;

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
     * @brief The type of the whole spline domain (cartesian product of 1D spline domain
     * and batch domain) preserving the underlying memory layout (order of dimensions).
     *
     * @tparam The batched discrete domain on which the interpolation points are defined.
     *
     * Example: For batched_interpolation_domain_type = DiscreteDomain<X,Y,Z> and a dimension of interest Y
     * (associated to a B-splines tag BSplinesY), this is DiscreteDomain<X,BSplinesY,Z>.
     */
    template <
            class BatchedInterpolationGrid,
            class = std::enable_if_t<ddc::is_discrete_domain_v<BatchedInterpolationGrid>>>
    using batched_basis_domain_type = ddc::replace_dim_of_t<
            BatchedInterpolationGrid,
            interpolation_discrete_dimension_type,
            basis_domain_type>;


public:
    /// @brief The number of equations defining the boundary condition at the lower bound.
    static constexpr int s_nbc_xmin = 0;

    /// @brief The number of equations defining the boundary condition at the upper bound.
    static constexpr int s_nbc_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

    /**
     * @brief Compute a spline approximation of a function.
     *
     * Use the values of a function (defined on
     * SplineBuilder::batched_interpolation_domain) and the derivatives of the
     * function at the boundaries (in the case of BoundCond::HERMITE only, defined
     * on SplineBuilder::batched_derivs_xmin_domain and SplineBuilder::batched_derivs_xmax_domain)
     * to calculate a spline approximation of this function.
     *
     * The spline approximation is stored as a ChunkSpan of coefficients
     * associated with B-splines.
     *
     * @param[out] spline The coefficients of the spline computed by this SplineBuilder.
     * @param[in] vals The values of the function on the interpolation mesh.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     * (used only with BoundCond::HERMITE lower boundary condition).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     * (used only with BoundCond::HERMITE upper boundary condition).
     */
    template <class Layout, class BatchedInterpolationGrid>
    void operator()(
            Field<DataType,
                  batched_basis_domain_type<BatchedInterpolationGrid>,
                  memory_space,
                  Layout> coeffs,
            ConstField<DataType, BatchedInterpolationGrid, Layout, memory_space> vals,
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
                = ddc::discrete_domain<Basis>().break_point_domain();
        ddc::parallel_deepcopy(coeffs[bp_idx_range], vals);
        if constexpr (Basis::is_uniform()) {
            IdxRange<basis_domain_type> extended_domain(
                    ddc::discrete_domain<Basis>().full_domain().remove_first(bp_idx_range));
            IdxRange<basis_domain_type>
                    repeat_domain(bp_idx_range.front(), extended_domain.extents());
            ddc::parallel_deepcopy(coeffs[extended_domain], vals[repeat_domain]);
        }
    }
};
