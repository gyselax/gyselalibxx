#pragma once

#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include <ddc/ddc.hpp>

#include "sll/math_tools.hpp"
#include "sll/matrix.hpp"
#include "sll/spline_boundary_conditions.hpp"
#include "sll/view.hpp"

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
class SplineBuilder
{
    static_assert(
            (BSplines::is_periodic() && (BcXmin == BoundCond::PERIODIC)
             && (BcXmax == BoundCond::PERIODIC))
            || (!BSplines::is_periodic() && (BcXmin != BoundCond::PERIODIC)
                && (BcXmax != BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename interpolation_mesh_type::continuous_dimension_type;

public:
    using bsplines_type = BSplines;

    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_mesh_type>;

public:
    static constexpr bool s_odd = BSplines::degree() % 2;

    static constexpr int s_nbe_xmin = n_boundary_equations(BcXmin, BSplines::degree());

    static constexpr int s_nbe_xmax = n_boundary_equations(BcXmax, BSplines::degree());

    static constexpr int s_nbc_xmin = n_user_input(BcXmin, BSplines::degree());

    static constexpr int s_nbc_xmax = n_user_input(BcXmax, BSplines::degree());

private:
    interpolation_domain_type m_interpolation_domain;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

    int m_offset;

public:
    SplineBuilder(interpolation_domain_type const& interpolation_domain);

    SplineBuilder(SplineBuilder const& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(SplineBuilder const& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    void operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> spline,
            ddc::ChunkSpan<double const, interpolation_domain_type> vals,
            std::optional<CDSpan1D> const derivs_xmin = std::nullopt,
            std::optional<CDSpan1D> const derivs_xmax = std::nullopt) const;

    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return m_interpolation_domain;
    }

    ddc::DiscreteDomain<BSplines> spline_domain() const noexcept
    {
        return ddc::discrete_space<BSplines>().full_domain();
    }

private:
    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(int lower_block_size, int upper_block_size);

    void compute_interpolant_degree1(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> spline,
            ddc::ChunkSpan<double const, interpolation_domain_type> vals) const;

    void build_matrix_system();
};

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::SplineBuilder(
        interpolation_domain_type const& interpolation_domain)
    : m_interpolation_domain(interpolation_domain)
    , m_dx((ddc::discrete_space<BSplines>().rmax() - ddc::discrete_space<BSplines>().rmin())
           / ddc::discrete_space<BSplines>().ncells())
    , matrix(nullptr)
    , m_offset(0)
{
    if constexpr (bsplines_type::is_periodic()) {
        // Calculate offset so that the matrix is diagonally dominant
        std::array<double, bsplines_type::degree() + 1> values_ptr;
        DSpan1D values(values_ptr.data(), bsplines_type::degree() + 1);
        ddc::DiscreteElement<interpolation_mesh_type> start(interpolation_domain.front());
        auto jmin = ddc::discrete_space<BSplines>()
                            .eval_basis(values, ddc::coordinate(start + BSplines::degree()));
        if constexpr (bsplines_type::degree() % 2 == 0) {
            m_offset = jmin.uid() - start.uid() + bsplines_type::degree() / 2 - BSplines::degree();
        } else {
            int const mid = bsplines_type::degree() / 2;
            m_offset = jmin.uid() - start.uid() + (values(mid) > values(mid + 1) ? mid : mid + 1)
                       - BSplines::degree();
        }
    }

    // Calculate block sizes
    int lower_block_size, upper_block_size;
    if constexpr (bsplines_type::is_uniform()) {
        compute_block_sizes_uniform(lower_block_size, upper_block_size);
    } else {
        compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
    }
    allocate_matrix(lower_block_size, upper_block_size);
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::compute_interpolant_degree1(
        ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> const spline,
        ddc::ChunkSpan<double const, interpolation_domain_type> const vals) const
{
    for (std::size_t i = 0; i < ddc::discrete_space<BSplines>().nbasis(); ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(i))
                = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
    }
    if constexpr (bsplines_type::is_periodic()) {
        spline(ddc::DiscreteElement<bsplines_type>(ddc::discrete_space<BSplines>().nbasis()))
                = spline(ddc::DiscreteElement<bsplines_type>(0));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::operator()(
        ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> const spline,
        ddc::ChunkSpan<double const, interpolation_domain_type> const vals,
        std::optional<CDSpan1D> const derivs_xmin,
        std::optional<CDSpan1D> const derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == ddc::discrete_space<BSplines>().nbasis() - s_nbe_xmin - s_nbe_xmax);
    // assert(spline.belongs_to_space(ddc::discrete_space<BSplines>()));
    // TODO: LOG Errors
    if constexpr (bsplines_type::degree() == 1)
        return compute_interpolant_degree1(spline, vals);

    assert((BcXmin == BoundCond::HERMITE)
           != (!derivs_xmin.has_value() || derivs_xmin->extent(0) == 0));
    assert((BcXmax == BoundCond::HERMITE)
           != (!derivs_xmax.has_value() || derivs_xmax->extent(0) == 0));

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmin == BoundCond::HERMITE) {
        assert(derivs_xmin->extent(0) == s_nbc_xmin);
        for (int i = s_nbc_xmin; i > 0; --i) {
            spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin - i))
                    = (*derivs_xmin)(i - 1) * ipow(m_dx, i + s_odd - 1);
        }
    }
    for (int i = s_nbc_xmin; i < s_nbc_xmin + m_offset; ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(i)) = 0.0;
    }

    for (int i = 0; i < m_interpolation_domain.extents(); ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin + i + m_offset))
                = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        assert(derivs_xmax->extent(0) == s_nbc_xmax);
        for (int i = 0; i < s_nbc_xmax; ++i) {
            spline(ddc::DiscreteElement<bsplines_type>(
                    ddc::discrete_space<BSplines>().nbasis() - s_nbc_xmax + i))
                    = (*derivs_xmax)(i)*ipow(m_dx, i + s_odd);
        }
    }

    DSpan1D const bcoef_section(spline.data() + m_offset, ddc::discrete_space<BSplines>().nbasis());
    matrix->solve_inplace(bcoef_section);

    if constexpr (bsplines_type::is_periodic()) {
        if (m_offset != 0) {
            for (int i = 0; i < m_offset; ++i) {
                spline(ddc::DiscreteElement<bsplines_type>(i))
                        = spline(ddc::DiscreteElement<bsplines_type>(
                                ddc::discrete_space<BSplines>().nbasis() + i));
            }
            for (std::size_t i = m_offset; i < bsplines_type::degree(); ++i) {
                spline(ddc::DiscreteElement<bsplines_type>(
                        ddc::discrete_space<BSplines>().nbasis() + i))
                        = spline(ddc::DiscreteElement<bsplines_type>(i));
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Compute num diags functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::compute_block_sizes_uniform(
        int& lower_block_size,
        int& upper_block_size) const
{
    switch (BcXmin) {
    case BoundCond::PERIODIC:
        upper_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::NATURAL:
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::NATURAL:
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::
        compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const
{
    switch (BcXmin) {
    case BoundCond::PERIODIC:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin + 1;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax + 1;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Initialize matrix functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::allocate_matrix(
        int lower_block_size,
        int upper_block_size)
{
    // Special case: linear spline
    // No need for matrix assembly
    if constexpr (bsplines_type::degree() == 1)
        return;

    int upper_band_width;
    if (bsplines_type::is_uniform()) {
        upper_band_width = bsplines_type::degree() / 2;
    } else {
        upper_band_width = bsplines_type::degree() - 1;
    }

    if constexpr (bsplines_type::is_periodic()) {
        matrix = Matrix::make_new_periodic_banded(
                ddc::discrete_space<BSplines>().nbasis(),
                upper_band_width,
                upper_band_width,
                bsplines_type::is_uniform());
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                ddc::discrete_space<BSplines>().nbasis(),
                upper_band_width,
                upper_band_width,
                bsplines_type::is_uniform(),
                upper_block_size,
                lower_block_size);
    }

    build_matrix_system();

    matrix->factorize();
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::build_matrix_system()
{
    // Hermite boundary conditions at xmin, if any
    if constexpr (BcXmin == BoundCond::HERMITE) {
        double derivs_ptr[(bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)];
        DSpan2D derivs(derivs_ptr, bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1);
        ddc::discrete_space<BSplines>()
                .eval_basis_and_n_derivs(derivs, discrete_space<BSplines>().rmin(), s_nbc_xmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        // iterate only to deg as last bspline is 0
        for (std::size_t i = 0; i < s_nbc_xmin; ++i) {
            for (std::size_t j = 0; j < bsplines_type::degree(); ++j) {
                matrix->set_element(i, j, derivs(j, s_nbc_xmin - i - 1 + s_odd));
            }
        }
    }

    // Interpolation points
    std::array<double, bsplines_type::degree() + 1> values_ptr;
    std::experimental::mdspan<double, std::experimental::extents<bsplines_type::degree() + 1>> const
            values(values_ptr.data());
    int start = m_interpolation_domain.front().uid();
    ddc::for_each(m_interpolation_domain, [&](auto ix) {
        auto jmin = ddc::discrete_space<BSplines>().eval_basis(
                values,
                coordinate(ddc::DiscreteElement<interpolation_mesh_type>(ix)));
        for (std::size_t s = 0; s < bsplines_type::degree() + 1; ++s) {
            int const j = modulo(
                    int(jmin.uid() - m_offset + s),
                    (int)ddc::discrete_space<BSplines>().nbasis());
            matrix->set_element(ix.uid() - start + s_nbc_xmin, j, values(s));
        }
    });

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == BoundCond::HERMITE) {
        std::array<double, (bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)>
                derivs_ptr;
        std::experimental::mdspan<
                double,
                std::experimental::
                        extents<bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1>> const
                derivs(derivs_ptr.data());

        ddc::discrete_space<BSplines>()
                .eval_basis_and_n_derivs(derivs, discrete_space<BSplines>().rmax(), s_nbc_xmax);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        int const i0 = ddc::discrete_space<BSplines>().nbasis() - s_nbc_xmax;
        int const j0 = ddc::discrete_space<BSplines>().nbasis() - bsplines_type::degree();
        for (std::size_t j = 0; j < bsplines_type::degree(); ++j) {
            for (std::size_t i = 0; i < s_nbc_xmax; ++i) {
                matrix->set_element(i0 + i, j0 + j, derivs(j + 1, i + s_odd));
            }
        }
    }
}
