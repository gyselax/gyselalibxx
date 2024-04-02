// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include "spline_builder.hpp"

namespace ddc {

/**
 * @brief A class for creating a 2D spline approximation of a function.
 *
 * A class which contains an operator () which can be used to build a 2D spline approximation
 * of a function. A 2D spline approximation uses a cross-product between two 1D spline builder.
 *
 * @see SplineBuilder
 */
template <
        class ExecSpace,
        class MemorySpace,
        class BSpline1,
        class BSpline2,
        class IDimI1,
        class IDimI2,
        ddc::BoundCond BcXmin1,
        ddc::BoundCond BcXmax1,
        ddc::BoundCond BcXmin2,
        ddc::BoundCond BcXmax2,
        ddc::SplineSolver Solver,
        class... IDimX>
class SplineBuilder2D
{
public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using builder_type1 = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSpline1,
            IDimI1,
            BcXmin1,
            BcXmax1,
            Solver,
            IDimX...>;
    using builder_type2 = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSpline2,
            IDimI2,
            BcXmin2,
            BcXmax2,
            Solver,
            std::conditional_t<std::is_same_v<IDimX, IDimI1>, BSpline1, IDimX>...>;
    using builder_deriv_type1 = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            BSpline1,
            IDimI1,
            BcXmin1,
            BcXmax1,
            Solver,
            std::conditional_t<
                    std::is_same_v<IDimX, IDimI2>,
                    typename builder_type2::deriv_type,
                    IDimX>...>;

private:
    /**
     * @brief Tag the dimension of the first 1D SplineBuilder.
     */
    using tag_type1 = typename builder_type1::bsplines_type::tag_type;
    /**
     * @brief Tag the dimension of the second 1D SplineBuilder.
     */
    using tag_type2 = typename builder_type2::bsplines_type::tag_type;

public:
    /**
     * @brief The type of the BSplines in the first dimension which are compatible with this class.
     */
    using bsplines_type1 = typename builder_type1::bsplines_type;
    /**
     * @brief The type of the BSplines in the second dimension which are compatible with this class.
     */
    using bsplines_type2 = typename builder_type2::bsplines_type;

    using deriv_type1 = typename builder_type1::deriv_type;
    using deriv_type2 = typename builder_type2::deriv_type;

    /**
     * @brief The type of the interpolation mesh in the first dimension used by this class.
     */
    using interpolation_mesh_type1 = typename builder_type1::interpolation_mesh_type;
    /**
     * @brief The type of the interpolation mesh in the second dimension used by this class.
     */
    using interpolation_mesh_type2 = typename builder_type2::interpolation_mesh_type;

    /**
     * @brief The type of the domain for the interpolation mesh is the first dimension used by this class.
     */
    using interpolation_domain_type1 = typename builder_type1::interpolation_mesh_type;
    /**
     * @brief The type of the domain for the interpolation mesh is the second dimension used by this class.
     */
    using interpolation_domain_type2 = typename builder_type2::interpolation_mesh_type;
    /**
     * @brief The type of the domain for the interpolation mesh is the 2D dimension used by this class.
     */
    using interpolation_domain_type
            = ddc::DiscreteDomain<interpolation_mesh_type1, interpolation_mesh_type2>;

    using vals_domain_type = ddc::DiscreteDomain<IDimX...>;

    using batch_domain_type
            = ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_remove_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type1, interpolation_mesh_type2>>>;

    using spline_domain_type
            = ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type1, interpolation_mesh_type2>,
                    ddc::detail::TypeSeq<bsplines_type1, bsplines_type2>>>;

    using derivs_domain_type1 = typename builder_type1::derivs_domain_type;
    using derivs_domain_type2
            = ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type2>,
                    ddc::detail::TypeSeq<deriv_type2>>>;
    using derivs_domain_type
            = ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type1, interpolation_mesh_type2>,
                    ddc::detail::TypeSeq<deriv_type1, deriv_type2>>>;

private:
    builder_type1 m_spline_builder1;
    builder_deriv_type1 m_spline_builder_deriv1;
    builder_type2 m_spline_builder2;

public:
    /**
     * @brief Create a new SplineBuilder2D.
     *
     * @param vals_domain
     *      The 2D domain on which points will be provided in order to
     *      create the 2D spline approximation.
     * @param cols_per_chunk The number of columns in the rhs passed to the underlying linear solver.
     * @param preconditionner_max_block_size The block size of in the block Jacobi preconditioner.
     */
    explicit SplineBuilder2D(
            vals_domain_type const& vals_domain,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : m_spline_builder1(vals_domain, cols_per_chunk, preconditionner_max_block_size)
        , m_spline_builder_deriv1(ddc::replace_dim_of<interpolation_mesh_type2, deriv_type2>(
                  m_spline_builder1.vals_domain(),
                  ddc::DiscreteDomain<deriv_type2>(
                          ddc::DiscreteElement<deriv_type2>(1),
                          ddc::DiscreteVector<deriv_type2>(bsplines_type2::degree() / 2))))
        , m_spline_builder2(
                  m_spline_builder1.spline_domain(),
                  cols_per_chunk,
                  preconditionner_max_block_size)
    {
    }

    /**
     * @brief Create a new SplineBuilder2D by copy
     *
     * @param x
     *      The SplineBuilder2D being copied.
     */
    SplineBuilder2D(SplineBuilder2D const& x) = delete;

    /**
     * @brief Create a new SplineBuilder2D by copy
     *
     * @param x
     *      The temporary SplineBuilder2D being copied.
     */
    SplineBuilder2D(SplineBuilder2D&& x) = default;

    ~SplineBuilder2D() = default;

    SplineBuilder2D& operator=(SplineBuilder2D const& x) = delete;


    /**
     * @brief Copy a SplineBuilder2D.
     *
     * @param x
     *      The temporary SplineBuilder2D being copied.
     * @returns A reference to this object.
     */
    SplineBuilder2D& operator=(SplineBuilder2D&& x) = default;

    vals_domain_type vals_domain() const noexcept
    {
        return m_spline_builder1.vals_domain();
    }

    /**
     * @brief Get the 2D dimension domain from which the approximation is defined.
     *
     * Get the 2D dimension  domain on which values of the function must be provided in order
     * to build a spline approximation of the function.
     *
     * @return The 2D dimension domain for the grid points.
     */
    interpolation_domain_type interpolation_domain() const noexcept
    {
        return ddc::DiscreteDomain<interpolation_domain_type1, interpolation_domain_type2>(
                m_spline_builder1.interpolation_domain(),
                m_spline_builder2.interpolation_domain());
    }

    batch_domain_type batch_domain() const noexcept
    {
        return ddc::remove_dims_of(vals_domain(), interpolation_domain());
    }

    /**
     * @brief Get the 2D domain on which the approximation is defined.
     *
     * Get the 2D domain of the basis-splines for which the coefficients of the spline
     * approximation must be calculated.
     *
     * @return The 2D domain for the splines.
     */
    ddc::DiscreteDomain<bsplines_type1, bsplines_type2> bsplines_domain()
            const noexcept // TODO : clarify name
    {
        return ddc::DiscreteDomain<bsplines_type1, bsplines_type2>(
                ddc::discrete_space<bsplines_type1>().full_domain(),
                ddc::discrete_space<bsplines_type2>().full_domain());
    }

    spline_domain_type spline_domain() const noexcept
    {
        return ddc::replace_dim_of<interpolation_mesh_type1, bsplines_type1>(
                ddc::replace_dim_of<
                        interpolation_mesh_type2,
                        bsplines_type2>(vals_domain(), bsplines_domain()),
                bsplines_domain());
    }

    /**
     * @brief Build a 2D spline approximation of a function.
     *
     * Use the values of a function at known grid points (as specified by
     * SplineBuilder2D::interpolation_domain_type) and the derivatives of the
     * function at the boundaries (if necessary for the chosen boundary
     * conditions) to calculate a 2D spline approximation of a function.
     *
     * The spline approximation is stored as a ChunkSpan of coefficients
     * associated with basis-splines.
     *
     * @param[out] spline
     *      The coefficients of the spline calculated by the function.
     * @param[in] vals
     *      The values of the function at the grid points.
     * @param[in] derivs_min1
     *      The values of the derivatives at the lower boundary in the first dimension.
     * @param[in] derivs_max1
     *      The values of the derivatives at the upper boundary in the first dimension.
     * @param[in] derivs_min2
     *      The values of the derivatives at the lower boundary in the second dimension.
     * @param[in] derivs_max2
     *      The values of the derivatives at the upper boundary in the second dimension.
     * @param[in] mixed_derivs_min1_min2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_max1_min2
     *      The values of the the cross-derivatives at the upper boundary in the first dimension
     *      and the lower boundary in the second dimension.
     * @param[in] mixed_derivs_min1_max2
     *      The values of the the cross-derivatives at the lower boundary in the first dimension
     *      and the upper boundary in the second dimension.
     * @param[in] mixed_derivs_max1_max2
     *      The values of the the cross-derivatives at the upper boundary in the first dimension
     *      and the upper boundary in the second dimension.
     */
    template <class Layout>
    void operator()(
            ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
            ddc::ChunkSpan<double const, vals_domain_type, Layout, memory_space> vals,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type1, Layout, memory_space>> const
                    derivs_min1
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type1, Layout, memory_space>> const
                    derivs_max1
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type2, Layout, memory_space>> const
                    derivs_min2
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type2, Layout, memory_space>> const
                    derivs_max2
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    mixed_derivs_min1_min2
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    mixed_derivs_max1_min2
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    mixed_derivs_min1_max2
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    mixed_derivs_max1_max2
            = std::nullopt) const;
};


template <
        class ExecSpace,
        class MemorySpace,
        class BSpline1,
        class BSpline2,
        class IDimI1,
        class IDimI2,
        ddc::BoundCond BcXmin1,
        ddc::BoundCond BcXmax1,
        ddc::BoundCond BcXmin2,
        ddc::BoundCond BcXmax2,
        ddc::SplineSolver Solver,
        class... IDimX>
template <class Layout>
void SplineBuilder2D<
        ExecSpace,
        MemorySpace,
        BSpline1,
        BSpline2,
        IDimI1,
        IDimI2,
        BcXmin1,
        BcXmax1,
        BcXmin2,
        BcXmax2,
        Solver,
        IDimX...>::
operator()(
        ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
        ddc::ChunkSpan<double const, vals_domain_type, Layout, memory_space> vals,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type1, Layout, memory_space>> const
                derivs_min1,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type1, Layout, memory_space>> const
                derivs_max1,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type2, Layout, memory_space>> const
                derivs_min2,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type2, Layout, memory_space>> const
                derivs_max2,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                mixed_derivs_min1_min2,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                mixed_derivs_max1_min2,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                mixed_derivs_min1_max2,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                mixed_derivs_max1_max2) const
{
    // TODO: perform computations along dimension 1 on different streams ?
    // Spline1-transform derivs_min2 (to spline1_deriv_min)
    ddc::Chunk spline1_deriv_min_alloc(
            m_spline_builder_deriv1.spline_domain(),
            ddc::KokkosAllocator<double, MemorySpace>());
    auto spline1_deriv_min = spline1_deriv_min_alloc.span_view();
    auto spline1_deriv_min_opt = std::optional(spline1_deriv_min.span_cview());
    if constexpr (BcXmin1 == ddc::BoundCond::HERMITE) {
        m_spline_builder_deriv1(
                spline1_deriv_min,
                *derivs_min2,
                mixed_derivs_min1_min2,
                mixed_derivs_max1_min2);
    } else {
        spline1_deriv_min_opt = std::nullopt;
    }

    // Spline1-transform vals (to spline1)
    ddc::Chunk spline1_alloc(
            m_spline_builder1.spline_domain(),
            ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline1 = spline1_alloc.span_view();

    m_spline_builder1(spline1, vals, derivs_min1, derivs_max1);

    // Spline1-transform derivs_max2 (to spline1_deriv_max)
    ddc::Chunk spline1_deriv_max_alloc(
            m_spline_builder_deriv1.spline_domain(),
            ddc::KokkosAllocator<double, MemorySpace>());
    auto spline1_deriv_max = spline1_deriv_max_alloc.span_view();
    auto spline1_deriv_max_opt = std::optional(spline1_deriv_max.span_cview());
    if constexpr (BcXmax1 == ddc::BoundCond::HERMITE) {
        m_spline_builder_deriv1(
                spline1_deriv_max,
                *derivs_max2,
                mixed_derivs_min1_max2,
                mixed_derivs_max1_max2);
    } else {
        spline1_deriv_max_opt = std::nullopt;
    }

    // Spline2-transform spline1
    m_spline_builder2(spline, spline1.span_cview(), spline1_deriv_min_opt, spline1_deriv_max_opt);
}
} // namespace ddc
