#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <ddc/Chunk>
#include <ddc/ChunkSpan>
#include <ddc/NonUniformDiscretization>
#include <ddc/UniformDiscretization>
#include <ddc/discretization>

#include "sll/math_tools.hpp"
#include "sll/matrix.hpp"

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
        std::exit(1);
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
            UniformDiscretization<tag_type>,
            NonUniformDiscretization<tag_type>>;

    using interpolation_domain_type = DiscreteDomain<interpolation_mesh_type>;

private:
    static constexpr bool s_odd = BSplines::degree() % 2;

    static constexpr int s_offset = BSplines::is_periodic() ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmin = BcXmin == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmax = BcXmin == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

private:
    std::unique_ptr<interpolation_domain_type> m_interpolation_domain;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

public:
    SplineBuilder();

    SplineBuilder(const SplineBuilder& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(const SplineBuilder& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    void operator()(
            ChunkSpan<double, DiscreteDomain<bsplines_type>> const& spline,
            ChunkSpan<double const, interpolation_domain_type> const& vals,
            DSpan1D const* derivs_xmin = nullptr,
            DSpan1D const* derivs_xmax = nullptr) const;

    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return *m_interpolation_domain;
    }

    DiscreteDomain<BSplines> spline_domain() const noexcept
    {
        return DiscreteDomain<BSplines>(
                DiscreteVector<BSplines>(discretization<BSplines>().size()));
    }

private:
    void compute_interpolation_points_uniform();

    void compute_interpolation_points_non_uniform();

    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(int kl, int ku);

    void compute_interpolant_degree1(
            ChunkSpan<double, bsplines_type>& spline,
            ChunkSpan<double, interpolation_domain_type> const& vals) const;

    void build_matrix_system();
};

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
SplineBuilder<BSplines, BcXmin, BcXmax>::SplineBuilder()
    : m_interpolation_domain(nullptr)
    , m_dx((discretization<BSplines>().rmax() - discretization<BSplines>().rmin())
           / discretization<BSplines>().ncells())
    , matrix(nullptr)
{
    int lower_block_size, upper_block_size;
    if constexpr (bsplines_type::is_uniform()) {
        compute_interpolation_points_uniform();
        compute_block_sizes_uniform(lower_block_size, upper_block_size);
    } else {
        compute_interpolation_points_non_uniform();
        compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
    }
    allocate_matrix(lower_block_size, upper_block_size);
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolant_degree1(
        ChunkSpan<double, bsplines_type>& spline,
        ChunkSpan<double, interpolation_domain_type> const& vals) const
{
    for (int i = 0; i < discretization<BSplines>().nbasis(); ++i) {
        spline(i) = vals(i);
    }
    if constexpr (bsplines_type::is_periodic()) {
        spline(discretization<BSplines>().nbasis()) = spline(0);
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::operator()(
        ChunkSpan<double, DiscreteDomain<bsplines_type>> const& spline,
        ChunkSpan<double const, interpolation_domain_type> const& vals,
        DSpan1D const* derivs_xmin,
        DSpan1D const* derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == discretization<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax);
    // assert(spline.belongs_to_space(discretization<BSplines>()));
    // TODO: LOG Errors
    if constexpr (bsplines_type::degree() == 1)
        return compute_interpolant_degree1(spline, vals);

    assert((BcXmin == BoundCond::HERMITE)
           != (derivs_xmin == nullptr || derivs_xmin->extent(0) == 0));
    assert((BcXmax == BoundCond::HERMITE)
           != (derivs_xmax == nullptr || derivs_xmax->extent(0) == 0));

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmin == BoundCond::HERMITE) {
        for (int i = s_nbc_xmin; i > 0; --i) {
            spline(DiscreteCoordinate<bsplines_type>(s_nbc_xmin - i))
                    = (*derivs_xmin)(i - 1) * ipow(m_dx, i + s_odd - 1);
        }
    }
    for (int i = s_nbc_xmin; i < s_nbc_xmin + s_offset; ++i) {
        spline(DiscreteCoordinate<bsplines_type>(i)) = 0.0;
    }

    for (int i = 0; i < m_interpolation_domain->extents(); ++i) {
        spline(DiscreteCoordinate<bsplines_type>(s_nbc_xmin + i + s_offset))
                = vals(DiscreteCoordinate<interpolation_mesh_type>(i));
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        for (int i = 0; i < s_nbc_xmax; ++i) {
            spline(DiscreteCoordinate<bsplines_type>(
                    discretization<BSplines>().nbasis() - s_nbc_xmax + i))
                    = (*derivs_xmax)(i)*ipow(m_dx, i + s_odd);
        }
    }

    DSpan1D bcoef_section(spline.data() + s_offset, discretization<BSplines>().nbasis());
    matrix->solve_inplace(bcoef_section);

    if constexpr (BcXmin == BoundCond::PERIODIC && s_offset != 0) {
        for (int i = 0; i < s_offset; ++i) {
            spline(DiscreteCoordinate<bsplines_type>(i)) = spline(
                    DiscreteCoordinate<bsplines_type>(discretization<BSplines>().nbasis() + i));
        }
        for (int i = s_offset; i < bsplines_type::degree(); ++i) {
            spline(DiscreteCoordinate<bsplines_type>(discretization<BSplines>().nbasis() + i))
                    = spline(DiscreteCoordinate<bsplines_type>(i));
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
    int const n_interp_pts = discretization<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax;

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        double constexpr shift = !s_odd ? 0.5 : 0.0;
        init_discretization<interpolation_mesh_type>(
                Coordinate<tag_type>(discretization<BSplines>().rmin() + shift * m_dx),
                Coordinate<tag_type>(m_dx));
        m_interpolation_domain = std::make_unique<interpolation_domain_type>(
                DiscreteVector<interpolation_mesh_type>(n_interp_pts));
    } else {
        std::vector<double> interp_pts(n_interp_pts);

        int n_iknots = n_interp_pts + bsplines_type::degree() - 1;
        std::vector<int> iknots(n_iknots);
        int i = 0;

        // Additional knots near x=xmin
        int n_to_fill_min = bsplines_type::degree() - s_nbc_xmin - 1;
        for (; i < n_to_fill_min; ++i) {
            if constexpr (BcXmin == BoundCond::GREVILLE)
                iknots[i] = 0;
            if constexpr (BcXmin == BoundCond::HERMITE)
                iknots[i] = -n_to_fill_min + i;
        }

        // Knots inside the domain
        for (int j = 0; j < discretization<BSplines>().ncells() + 1; ++i, ++j) {
            iknots[i] = j;
        }

        // Additional knots near x=xmax
        for (int j = 1; i < n_iknots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                iknots[i] = discretization<BSplines>().ncells();
            if constexpr (BcXmax == BoundCond::HERMITE)
                iknots[i] = discretization<BSplines>().ncells() + j;
        }

        for (int j = 0; j < n_interp_pts; ++j) {
            int isum = sum(iknots.data() + j, bsplines_type::degree());
            interp_pts[j]
                    = discretization<BSplines>().rmin() + m_dx * isum / bsplines_type::degree();
        }

        // Non-periodic case, odd degree: fix round-off issues
        if constexpr (s_odd) {
            interp_pts[0] = discretization<BSplines>().rmin();
            interp_pts[n_interp_pts - 1] = discretization<BSplines>().rmax();
        }
        init_discretization<interpolation_mesh_type>(interp_pts);
        m_interpolation_domain = std::make_unique<interpolation_domain_type>(
                DiscreteVector<NonUniformDiscretization<tag_type>>(interp_pts.size()));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolation_points_non_uniform()
{
    int n_interp_pts = discretization<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax;
    std::vector<double> interp_pts(n_interp_pts);

    int n_temp_knots = n_interp_pts - 1 + bsplines_type::degree();
    double temp_knots[n_temp_knots];

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        for (int i = 0; i < n_interp_pts - 1 + bsplines_type::degree(); ++i) {
            temp_knots[i] = discretization<BSplines>().get_knot(
                    1 - bsplines_type::degree() + s_offset + i);
        }
    } else {
        int i = 0;
        int n_start_pts = bsplines_type::degree() - s_nbc_xmin - 1;

        // Initialise knots relevant to the xmin boundary condition
        for (; i < n_start_pts; ++i) {
            // As xmin_bc is a const variable the compiler should optimize
            // for(if..else..) to if(for..)else(for...)
            if constexpr (BcXmin == BoundCond::GREVILLE)
                temp_knots[i] = discretization<BSplines>().get_knot(0);
            if constexpr (BcXmin == BoundCond::HERMITE)
                temp_knots[i] = 2.0 * discretization<BSplines>().get_knot(0)
                                - discretization<BSplines>().get_knot(n_start_pts - i);
        }

        // Initialise central knots
        for (int j = 0; j < discretization<BSplines>().npoints(); ++i, ++j) {
            temp_knots[i] = discretization<BSplines>().get_knot(j);
        }

        // Initialise knots relevant to the xmax boundary condition
        for (int j = 0; i < n_temp_knots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                temp_knots[i]
                        = discretization<BSplines>().get_knot(discretization<BSplines>().ncells());
            if constexpr (BcXmax == BoundCond::HERMITE)
                temp_knots[i] = 2.0
                                        * discretization<BSplines>().get_knot(
                                                discretization<BSplines>().ncells())
                                - discretization<BSplines>().get_knot(
                                        discretization<BSplines>().ncells() - 1 - j);
        }
    }

    // Compute interpolation points using Greville-style averaging
    double inv_deg = 1.0 / bsplines_type::degree();
    for (int i = 0; i < n_interp_pts; ++i) {
        interp_pts[i] = sum(temp_knots + i, bsplines_type::degree()) * inv_deg;
    }

    // Periodic case: apply periodic BCs to interpolation points
    if constexpr (BcXmin == BoundCond::PERIODIC) {
        double zone_width(discretization<BSplines>().rmax() - discretization<BSplines>().rmin());
        for (int i = 0; i < n_interp_pts; ++i) {
            interp_pts[i] = modulo(interp_pts[i] - s_nbc_xmin, zone_width)
                            + discretization<BSplines>().rmin();
        }
    }
    // Non-periodic case, odd degree: fix round-off issues
    else {
        if constexpr (s_odd) {
            interp_pts[0] = discretization<BSplines>().rmin();
            interp_pts[n_interp_pts - 1] = discretization<BSplines>().rmax();
        }
    }

    init_discretization<interpolation_mesh_type>(interp_pts);
    m_interpolation_domain = std::make_unique<interpolation_domain_type>(
            DiscreteCoordinate<NonUniformDiscretization<tag_type>>(interp_pts.size()));
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
        upper_block_size = (bsplines_type::degree() + 1) / 2;
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
        lower_block_size = (bsplines_type::degree() + 1) / 2;
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
        upper_band_width = (bsplines_type::degree() + 1) / 2;
    }

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        matrix = Matrix::make_new_periodic_banded(
                discretization<BSplines>().nbasis(),
                upper_band_width,
                upper_band_width,
                bsplines_type::is_uniform());
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                discretization<BSplines>().nbasis(),
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
        discretization<BSplines>().eval_basis_and_n_derivs(
                discretization<BSplines>().rmin(),
                s_nbc_xmin,
                derivs,
                jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (int j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        // iterate only to deg as last bspline is 0
        for (int j = 0; j < s_nbc_xmin; ++j) {
            for (int i = 0; i < bsplines_type::degree(); ++i) {
                // Elements are set in Fortran order as they are LAPACK input
                matrix->set_element(i, j, derivs(i, s_nbc_xmin - j - 1 + s_odd));
            }
        }
    }

    // Interpolation points
    double values_ptr[bsplines_type::degree() + 1];
    DSpan1D values(values_ptr, bsplines_type::degree() + 1);
    for (int i = 0; i < discretization<BSplines>().nbasis() - s_nbc_xmin - s_nbc_xmax; ++i) {
        discretization<BSplines>()
                .eval_basis(to_real(DiscreteCoordinate<interpolation_mesh_type>(i)), values, jmin);
        for (int s = 0; s < bsplines_type::degree() + 1; ++s) {
            int j = modulo(jmin - s_offset + s, (int)discretization<BSplines>().nbasis());
            matrix->set_element(j, i + s_nbc_xmin, values(s));
        }
    }

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == BoundCond::HERMITE) {
        double derivs_ptr[(bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)];
        DSpan2D derivs(derivs_ptr, bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1);

        discretization<BSplines>().eval_basis_and_n_derivs(
                discretization<BSplines>().rmax(),
                s_nbc_xmax,
                derivs,
                jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (int j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        int i0 = discretization<BSplines>().nbasis() - bsplines_type::degree();
        int j0 = discretization<BSplines>().nbasis() - s_nbc_xmax;
        for (int i = 0; i < bsplines_type::degree(); ++i) {
            for (int j = 0; j < s_nbc_xmax; ++j) {
                matrix->set_element(i0 + i, j0 + j, derivs(i + 1, j + s_odd));
            }
        }
    }
}
