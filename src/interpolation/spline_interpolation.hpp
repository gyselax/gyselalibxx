// SPDX-License-Identifier: MIT

#pragma once

#include <optional>
#include <type_traits>

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

/**
 * @brief A class for building spline coefficients that satisfies the IInterpolationBuilder interface.
 *
 * This class wraps ddc::SplineBuilder and exposes it through the IInterpolationBuilder interface.
 *
 * @tparam ExecSpace The Kokkos execution space on which the spline approximation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data is stored.
 * @tparam DataType The data type (must be double, matching ddc::SplineBuilder).
 * @tparam BSplines The discrete dimension representing the B-splines.
 * @tparam InterpolationGrid The discrete dimension on which interpolation points are defined.
 * @tparam BcLower The lower boundary condition.
 * @tparam BcUpper The upper boundary condition.
 * @tparam Solver The SplineSolver giving the backend used to perform the spline approximation.
 * @tparam BatchedInterpolationIdxRange The batched discrete domain on which interpolation is defined.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class BSplines,
        class InterpolationGrid,
        ddc::BoundCond BcLower,
        ddc::BoundCond BcUpper,
        ddc::SplineSolver Solver,
        class BatchedInterpolationIdxRange = IdxRange<InterpolationGrid>>
class SplineBuilder1D
    : public IInterpolationBuilder<
              ExecSpace,
              MemorySpace,
              DataType,
              InterpolationGrid,
              BSplines,
              BatchedInterpolationIdxRange>
{
    static_assert(
            std::is_same_v<DataType, double>,
            "SplineBuilder1D only supports double as DataType, matching ddc::SplineBuilder.");

    using base_type = IInterpolationBuilder<
            ExecSpace,
            MemorySpace,
            DataType,
            InterpolationGrid,
            BSplines,
            BatchedInterpolationIdxRange>;

    using builder_type = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSplines,
            InterpolationGrid,
            BcLower,
            BcUpper,
            Solver>;

    builder_type const& m_builder;

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using typename base_type::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using typename base_type::memory_space;

    /// @brief The type of the interpolation continuous dimension used by this class.
    using typename base_type::continuous_dimension_type;

    /// @brief The type of the interpolation discrete dimension used by this class.
    using typename base_type::interpolation_grid_type;

    /// @brief The type of the domain for the 1D interpolation mesh used by this class.
    using typename base_type::interpolation_idx_range_type;

    /// @brief The grid on which the interpolation coefficients should be provided (equals BSplines).
    using typename base_type::basis_domain_type;

    /// @brief The discrete dimension representing the B-splines (same as basis_domain_type).
    using bsplines_type = BSplines;

    /// @brief The number of equations defining the boundary condition at the lower bound.
    static constexpr int s_nbc_xmin = builder_type::s_nbc_xmin;

    /// @brief The number of equations defining the boundary condition at the upper bound.
    static constexpr int s_nbc_xmax = builder_type::s_nbc_xmax;

    /// @brief The boundary condition implemented at the lower bound.
    static constexpr ddc::BoundCond s_bc_xmin = BcLower;

    /// @brief The boundary condition implemented at the upper bound.
    static constexpr ddc::BoundCond s_bc_xmax = BcUpper;

    /// @brief The SplineSolver giving the backend used to perform the spline approximation.
    static constexpr ddc::SplineSolver s_spline_solver = Solver;

public:
    /**
     * @brief Construct a SplineBuilder1D acting on the given 1D interpolation domain.
     *
     * @param interpolation_domain The 1D domain on which the interpolation points are defined.
     * @param cols_per_chunk Optional solver parameter (see ddc::SplineBuilder).
     * @param preconditioner_max_block_size Optional solver parameter (see ddc::SplineBuilder).
     */
    explicit SplineBuilder1D(builder_type const& builder) : m_builder(builder) {}

    /// @brief Copy-constructor is deleted.
    SplineBuilder1D(SplineBuilder1D const&) = delete;

    /// @brief Move-constructor.
    SplineBuilder1D(SplineBuilder1D&&) = default;

    /// @brief Destructor.
    ~SplineBuilder1D() = default;

    /// @brief Copy-assignment is deleted.
    SplineBuilder1D& operator=(SplineBuilder1D const&) = delete;

    /// @brief Move-assignment.
    SplineBuilder1D& operator=(SplineBuilder1D&&) = default;

    /**
     * @brief Get the 1D domain for the interpolation mesh.
     *
     * @return The 1D interpolation mesh domain.
     */
    interpolation_idx_range_type interpolation_domain() const noexcept
    {
        return m_builder.interpolation_domain();
    }

    /**
     * @brief Compute the spline coefficients for a function.
     *
     * Delegates to ddc::SplineBuilder to compute B-spline coefficients from function values.
     *
     * @param[out] coeffs The B-spline coefficients computed by this builder.
     * @param[in] vals The values of the function at the interpolation points.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     *                  (used only with BoundCond::HERMITE lower boundary condition).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     *                  (used only with BoundCond::HERMITE upper boundary condition).
     */
    void operator()(
            Field<DataType, typename base_type::batched_basis_idx_range_type, memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<
                    DataType,
                    typename base_type::batched_derivs_idx_range_type,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    typename base_type::batched_derivs_idx_range_type,
                    memory_space>> derivs_xmax
            = std::nullopt) const final
    {
        m_builder(coeffs, vals);
    }

    /**
     * @brief Get the whole domain on which derivatives on lower boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    typename base_type::batched_derivs_idx_range_type batched_derivs_xmin_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept final
    {
        return m_builder.batched_derivs_xmin_domain(batched_interpolation_domain);
    }

    /**
     * @brief Get the whole domain on which derivatives on upper boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    typename base_type::batched_derivs_idx_range_type batched_derivs_xmax_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept final
    {
        return m_builder.batched_derivs_xmax_domain(batched_interpolation_domain);
    }
};
