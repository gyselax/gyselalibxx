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
#include "sll/view.hpp"

enum class BoundCond {
    // Periodic boundary condition u(1)=u(n)
    PERIODIC,
    // Hermite boundary condition
    HERMITE,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    GREVILLE,
};

static inline std::ostream& operator<<(std::ostream& out, BoundCond bc)
{
    switch (bc) {
    case BoundCond::PERIODIC:
        return out << "PERIODIC";
    case BoundCond::HERMITE:
        return out << "HERMITE";
    case BoundCond::GREVILLE:
        return out << "GREVILLE";
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
class SplineBuilder
{
    static_assert(
            (BSplines::is_periodic() && (BcXmin == BoundCond::PERIODIC)
             && (BcXmax == BoundCond::PERIODIC))
            || (!BSplines::is_periodic() && (BcXmin != BoundCond::PERIODIC)
                && (BcXmax != BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename BSplines::tag_type;

public:
    using bsplines_type = BSplines;

    // No need to check boundary conditions, it shall fail if it periodic with non-periodic boundary conditions
    using interpolation_mesh_type = std::conditional_t<
            BSplines::is_uniform() && BSplines::is_periodic(),
            UniformPointSampling<tag_type>,
            NonUniformPointSampling<tag_type>>;

    using interpolation_domain_type = DiscreteDomain<interpolation_mesh_type>;

private:
    static constexpr bool s_odd = BSplines::degree() % 2;

    static constexpr int s_offset = BSplines::is_periodic() ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmin = BcXmin == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmax = BcXmax == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

private:
    std::unique_ptr<interpolation_domain_type> m_interpolation_domain;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

public:
    SplineBuilder();

    SplineBuilder(SplineBuilder const& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(SplineBuilder const& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    void operator()(
            ChunkSpan<double, DiscreteDomain<bsplines_type>> spline,
            ChunkSpan<double const, interpolation_domain_type> vals,
            std::optional<DSpan1D> const derivs_xmin = std::nullopt,
            std::optional<DSpan1D> const derivs_xmax = std::nullopt) const;

    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return *m_interpolation_domain;
    }

    DiscreteDomain<BSplines> spline_domain() const noexcept
    {
        return discrete_space<BSplines>().full_domain();
    }

private:
    void compute_interpolation_points_uniform();

    int compute_interpolation_points_non_uniform();

    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(int lower_block_size, int upper_block_size, int diag_shift);

    void compute_interpolant_degree1(
            ChunkSpan<double, DiscreteDomain<bsplines_type>> spline,
            ChunkSpan<double const, interpolation_domain_type> vals) const;

    void build_matrix_system();
};

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
SplineBuilder<BSplines, BcXmin, BcXmax>::SplineBuilder()
    : m_interpolation_domain(nullptr)
    , m_dx((discrete_space<BSplines>().rmax() - discrete_space<BSplines>().rmin())
           / discrete_space<BSplines>().ncells())
    , matrix(nullptr)
{
    int diag_shift = 0;
    int lower_block_size, upper_block_size;
    if constexpr (bsplines_type::is_uniform()) {
        compute_interpolation_points_uniform();
        compute_block_sizes_uniform(lower_block_size, upper_block_size);
    } else {
        diag_shift = compute_interpolation_points_non_uniform();
        compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
    }
    allocate_matrix(lower_block_size, upper_block_size, diag_shift);
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolant_degree1(
        ChunkSpan<double, DiscreteDomain<bsplines_type>> const spline,
        ChunkSpan<double const, interpolation_domain_type> const vals) const
{
    for (std::size_t i = 0; i < discrete_space<BSplines>().nbasis(); ++i) {
        spline(DiscreteElement<bsplines_type>(i))
                = vals(DiscreteElement<interpolation_mesh_type>(i));
    }
    if constexpr (bsplines_type::is_periodic()) {
        spline(DiscreteElement<bsplines_type>(discrete_space<BSplines>().nbasis()))
                = spline(DiscreteElement<bsplines_type>(0));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::operator()(
        ChunkSpan<double, DiscreteDomain<bsplines_type>> const spline,
        ChunkSpan<double const, interpolation_domain_type> const vals,
        std::optional<DSpan1D> const derivs_xmin,
        std::optional<DSpan1D> const derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == discrete_space<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax);
    // assert(spline.belongs_to_space(discrete_space<BSplines>()));
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
        for (int i = s_nbc_xmin; i > 0; --i) {
            spline(DiscreteElement<bsplines_type>(s_nbc_xmin - i))
                    = (*derivs_xmin)(i - 1) * ipow(m_dx, i + s_odd - 1);
        }
    }
    for (int i = s_nbc_xmin; i < s_nbc_xmin + s_offset; ++i) {
        spline(DiscreteElement<bsplines_type>(i)) = 0.0;
    }

    for (int i = 0; i < m_interpolation_domain->extents(); ++i) {
        spline(DiscreteElement<bsplines_type>(s_nbc_xmin + i + s_offset))
                = vals(DiscreteElement<interpolation_mesh_type>(i));
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        for (int i = 0; i < s_nbc_xmax; ++i) {
            spline(DiscreteElement<bsplines_type>(
                    discrete_space<BSplines>().nbasis() - s_nbc_xmax + i))
                    = (*derivs_xmax)(i)*ipow(m_dx, i + s_odd);
        }
    }

    DSpan1D const bcoef_section(spline.data() + s_offset, discrete_space<BSplines>().nbasis());
    matrix->solve_inplace(bcoef_section);

    if constexpr (bsplines_type::is_periodic() && s_offset != 0) {
        for (std::size_t i = 0; i < s_offset; ++i) {
            spline(DiscreteElement<bsplines_type>(i)) = spline(
                    DiscreteElement<bsplines_type>(discrete_space<BSplines>().nbasis() + i));
        }
        for (std::size_t i = s_offset; i < bsplines_type::degree(); ++i) {
            spline(DiscreteElement<bsplines_type>(discrete_space<BSplines>().nbasis() + i))
                    = spline(DiscreteElement<bsplines_type>(i));
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                  Compute interpolation points functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolation_points_uniform()
{
    std::size_t const n_interp_pts = discrete_space<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax;

    if constexpr (bsplines_type::is_periodic()) {
        double constexpr shift = !s_odd ? 0.5 : 0.0;
        init_discrete_space<interpolation_mesh_type>(
                Coordinate<tag_type>(discrete_space<BSplines>().rmin() + shift * m_dx),
                Coordinate<tag_type>(m_dx));
        m_interpolation_domain = std::make_unique<interpolation_domain_type>(
                DiscreteVector<interpolation_mesh_type>(n_interp_pts));
    } else {
        std::vector<double> interp_pts(n_interp_pts);

        std::size_t const n_iknots = n_interp_pts + bsplines_type::degree() - 1;
        std::vector<int> iknots(n_iknots);
        std::size_t i = 0;

        // Additional knots near x=xmin
        int const n_to_fill_min = bsplines_type::degree() - s_nbc_xmin - 1;
        for (; i < n_to_fill_min; ++i) {
            if constexpr (BcXmin == BoundCond::GREVILLE)
                iknots[i] = 0;
            if constexpr (BcXmin == BoundCond::HERMITE)
                iknots[i] = -n_to_fill_min + i;
        }

        // Knots inside the domain
        for (std::size_t j = 0; j < discrete_space<BSplines>().ncells() + 1; ++i, ++j) {
            iknots[i] = j;
        }

        // Additional knots near x=xmax
        for (std::size_t j = 1; i < n_iknots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                iknots[i] = discrete_space<BSplines>().ncells();
            if constexpr (BcXmax == BoundCond::HERMITE)
                iknots[i] = discrete_space<BSplines>().ncells() + j;
        }

        for (std::size_t j = 0; j < n_interp_pts; ++j) {
            int const isum = sum(iknots.data() + j, bsplines_type::degree());
            interp_pts[j]
                    = discrete_space<BSplines>().rmin() + m_dx * isum / bsplines_type::degree();
        }

        // Non-periodic case, odd degree: fix round-off issues
        if constexpr (s_odd) {
            interp_pts[0] = discrete_space<BSplines>().rmin();
            interp_pts[n_interp_pts - 1] = discrete_space<BSplines>().rmax();
        }
        init_discrete_space<interpolation_mesh_type>(interp_pts);
        m_interpolation_domain = std::make_unique<interpolation_domain_type>(
                DiscreteVector<NonUniformPointSampling<tag_type>>(interp_pts.size()));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
int SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolation_points_non_uniform()
{
    int diag_shift = 0;
    std::size_t const n_interp_pts = discrete_space<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax;
    std::vector<double> interp_pts(n_interp_pts);

    std::size_t const n_temp_knots = n_interp_pts - 1 + bsplines_type::degree();
    std::vector<double> temp_knots(n_temp_knots);

    if constexpr (bsplines_type::is_periodic()) {
        for (std::size_t i = 0; i < n_interp_pts - 1 + bsplines_type::degree(); ++i) {
            temp_knots[i] = discrete_space<BSplines>().get_knot(
                    1 - bsplines_type::degree() + s_offset + i);
        }
    } else {
        std::size_t i = 0;
        std::size_t const n_start_pts = bsplines_type::degree() - s_nbc_xmin - 1;

        // Initialise knots relevant to the xmin boundary condition
        for (; i < n_start_pts; ++i) {
            // As xmin_bc is a const variable the compiler should optimize
            // for(if..else..) to if(for..)else(for...)
            if constexpr (BcXmin == BoundCond::GREVILLE)
                temp_knots[i] = discrete_space<BSplines>().get_knot(0);
            if constexpr (BcXmin == BoundCond::HERMITE)
                temp_knots[i] = 2.0 * discrete_space<BSplines>().get_knot(0)
                                - discrete_space<BSplines>().get_knot(n_start_pts - i);
        }

        // Initialise central knots
        for (std::size_t j = 0; j < discrete_space<BSplines>().npoints(); ++i, ++j) {
            temp_knots[i] = discrete_space<BSplines>().get_knot(j);
        }

        // Initialise knots relevant to the xmax boundary condition
        for (std::size_t j = 0; i < n_temp_knots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                temp_knots[i]
                        = discrete_space<BSplines>().get_knot(discrete_space<BSplines>().ncells());
            if constexpr (BcXmax == BoundCond::HERMITE)
                temp_knots[i] = 2.0
                                        * discrete_space<BSplines>().get_knot(
                                                discrete_space<BSplines>().ncells())
                                - discrete_space<BSplines>().get_knot(
                                        discrete_space<BSplines>().ncells() - 1 - j);
        }
    }

    // Compute interpolation points using Greville-style averaging
    double const inv_deg = 1.0 / bsplines_type::degree();
    for (std::size_t i = 0; i < n_interp_pts; ++i) {
        interp_pts[i] = sum(temp_knots.data() + i, bsplines_type::degree()) * inv_deg;
    }

    // Periodic case: apply periodic BCs to interpolation points
    if constexpr (bsplines_type::is_periodic()) {
        double const zone_width
                = discrete_space<BSplines>().rmax() - discrete_space<BSplines>().rmin();
        if (interp_pts[0] < discrete_space<BSplines>().rmin()) {
            // Count the number of interpolation points that need shifting to preserve the ordering
            while (interp_pts[diag_shift] < discrete_space<BSplines>().rmin()) {
                temp_knots[diag_shift] = modulo(interp_pts[diag_shift] - s_nbc_xmin, zone_width)
                                         + discrete_space<BSplines>().rmin();
                diag_shift++;
            }
            // Shift the points
            for (std::size_t i = 0; i < n_interp_pts - diag_shift; ++i) {
                interp_pts[i] = interp_pts[i + diag_shift];
            }
            for (std::size_t i = 0; i < diag_shift; ++i) {
                interp_pts[n_interp_pts - diag_shift + i] = temp_knots[i];
            }
        } else if (interp_pts[n_interp_pts - 1] > discrete_space<BSplines>().rmax()) {
            // Count the number of interpolation points that need shifting to preserve the ordering
            while (interp_pts[n_interp_pts - 1 + diag_shift] > discrete_space<BSplines>().rmin()) {
                temp_knots[-diag_shift]
                        = modulo(interp_pts[n_interp_pts - 1 + diag_shift] - s_nbc_xmin, zone_width)
                          + discrete_space<BSplines>().rmin();
                diag_shift--;
            }
            // Shift the points
            for (std::size_t i = 0; i < n_interp_pts + diag_shift; ++i) {
                interp_pts[i - diag_shift] = interp_pts[i];
            }
            for (std::size_t i = 0; i < -diag_shift; ++i) {
                interp_pts[-diag_shift - 1 - i] = temp_knots[i];
            }
        }
    }
    // Non-periodic case, odd degree: fix round-off issues
    else {
        if constexpr (s_odd) {
            interp_pts[0] = discrete_space<BSplines>().rmin();
            interp_pts[n_interp_pts - 1] = discrete_space<BSplines>().rmax();
        }
    }

    init_discrete_space<interpolation_mesh_type>(interp_pts);
    m_interpolation_domain = std::make_unique<interpolation_domain_type>(
            DiscreteVector<NonUniformPointSampling<tag_type>>(interp_pts.size()));

    return diag_shift;
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Compute num diags functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_block_sizes_uniform(
        int& lower_block_size,
        int& upper_block_size) const
{
    switch (BcXmin) {
    case BoundCond::PERIODIC:
        upper_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        break; // TODO: throw error
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        break; // TODO: throw error
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_block_sizes_non_uniform(
        int& lower_block_size,
        int& upper_block_size) const
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
        break; // TODO: throw error
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
        break; // TODO: throw error
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Initialize matrix functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::allocate_matrix(
        int lower_block_size,
        int upper_block_size,
        int diag_shift)
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
                discrete_space<BSplines>().nbasis(),
                upper_band_width - diag_shift,
                upper_band_width + diag_shift,
                bsplines_type::is_uniform());
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                discrete_space<BSplines>().nbasis(),
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

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::build_matrix_system()
{
    int jmin;

    // Hermite boundary conditions at xmin, if any
    if constexpr (BcXmin == BoundCond::HERMITE) {
        double derivs_ptr[(bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)];
        DSpan2D derivs(derivs_ptr, bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1);
        discrete_space<BSplines>().eval_basis_and_n_derivs(
                derivs,
                jmin,
                discrete_space<BSplines>().rmin(),
                s_nbc_xmin);

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
    for (std::size_t i = 0; i < discrete_space<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax;
         ++i) {
        discrete_space<BSplines>()
                .eval_basis(values, jmin, coordinate(DiscreteElement<interpolation_mesh_type>(i)));
        for (std::size_t s = 0; s < bsplines_type::degree() + 1; ++s) {
            int const j
                    = modulo(int(jmin - s_offset + s), (int)discrete_space<BSplines>().nbasis());
            matrix->set_element(i + s_nbc_xmin, j, values(s));
        }
    }

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == BoundCond::HERMITE) {
        std::array<double, (bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)>
                derivs_ptr;
        std::experimental::mdspan<
                double,
                std::experimental::
                        extents<bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1>> const
                derivs(derivs_ptr.data());

        discrete_space<BSplines>().eval_basis_and_n_derivs(
                derivs,
                jmin,
                discrete_space<BSplines>().rmax(),
                s_nbc_xmax);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        int const i0 = discrete_space<BSplines>().nbasis() - s_nbc_xmax;
        int const j0 = discrete_space<BSplines>().nbasis() - bsplines_type::degree();
        for (std::size_t j = 0; j < bsplines_type::degree(); ++j) {
            for (std::size_t i = 0; i < s_nbc_xmax; ++i) {
                matrix->set_element(i0 + i, j0 + j, derivs(j + 1, i + s_odd));
            }
        }
    }
}
