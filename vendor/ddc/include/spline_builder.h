#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include <math_tools.h>
#include <matrix.h>

#include "block.h"
#include "blockview.h"
#include "blockview_spline.h"
#include "mdomain.h"
#include "non_uniform_mesh.h"
#include "uniform_mesh.h"

enum class BoundCond {
    // Periodic boundary condition u(1)=u(n)
    PERIODIC,
    // Hermite boundary condition
    HERMITE,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    GREVILLE,
};

static std::ostream& operator<<(std::ostream& out, BoundCond bc)
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
            BSplines::is_periodic()
            == ((BcXmin == BoundCond::PERIODIC) == (BcXmax == BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename BSplines::tag_type;

public:
    using bsplines_type = BSplines;

    // No need to check boundary conditions, it shall fail if it periodic with non-periodic boundary conditions
    using interpolation_domain_type = std::conditional_t<
            BSplines::is_uniform() && BSplines::is_periodic(),
            MDomainImpl<UniformMesh<tag_type>>,
            MDomainImpl<NonUniformMesh<tag_type>>>;

private:
    static constexpr bool s_odd = BSplines::degree() % 2;

    static constexpr int s_offset = BSplines::is_periodic() ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmin = BcXmin == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

    static constexpr int s_nbc_xmax = BcXmin == BoundCond::HERMITE ? BSplines::degree() / 2 : 0;

public:
    static int compute_num_cells(int degree, BoundCond xmin, BoundCond xmax, int nipts);

private:
    std::unique_ptr<interpolation_domain_type> m_interpolation_domain;

    bsplines_type const& m_bsplines;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

public:
    SplineBuilder() = delete;

    SplineBuilder(BSplines const& bsplines);

    SplineBuilder(const SplineBuilder& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(const SplineBuilder& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    void operator()(
            BlockView<bsplines_type, double>& spline,
            BlockView<interpolation_domain_type, double const> const& vals,
            DSpan1D const* derivs_xmin = nullptr,
            DSpan1D const* derivs_xmax = nullptr) const;

    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return *m_interpolation_domain;
    }

private:
    void compute_interpolation_points_uniform();

    void compute_interpolation_points_non_uniform();

    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(int kl, int ku);

    void compute_interpolant_degree1(
            BlockView<bsplines_type, double>& spline,
            BlockView<interpolation_domain_type, double> const& vals) const;

    void build_matrix_system();
};

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
SplineBuilder<BSplines, BcXmin, BcXmax>::SplineBuilder(BSplines const& bsplines)
    : m_interpolation_domain(nullptr)
    , m_bsplines(bsplines)
    , m_dx((bsplines.rmax() - bsplines.rmin()) / bsplines.ncells())
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
        BlockView<bsplines_type, double>& spline,
        BlockView<interpolation_domain_type, double> const& vals) const
{
    for (int i(0); i < m_bsplines.nbasis(); ++i) {
        spline(i) = vals(i);
    }
    if constexpr (bsplines_type::is_periodic()) {
        spline(m_bsplines.nbasis()) = spline(0);
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::operator()(
        BlockView<bsplines_type, double>& spline,
        BlockView<interpolation_domain_type, double const> const& vals,
        DSpan1D const* derivs_xmin,
        DSpan1D const* derivs_xmax) const
{
    assert(vals.extent(0) == m_bsplines.nbasis() - s_nbc_xmin - s_nbc_xmax);
    // assert(spline.belongs_to_space(m_bsplines));
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
        for (int i(s_nbc_xmin); i > 0; --i) {
            spline(s_nbc_xmin - i) = (*derivs_xmin)(i - 1) * ipow(m_dx, i + s_odd - 1);
        }
    }
    for (int i(s_nbc_xmin); i < s_nbc_xmin + s_offset; ++i) {
        spline(i) = 0.0;
    }

    for (int i(0); i < m_interpolation_domain->template extent<tag_type>(); ++i) {
        spline(s_nbc_xmin + i + s_offset) = vals(i);
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        for (int i(0); i < s_nbc_xmax; ++i) {
            spline(m_bsplines.nbasis() - s_nbc_xmax + i) = (*derivs_xmax)(i)*ipow(m_dx, i + s_odd);
        }
    }

    DSpan1D bcoef_section(spline.raw_view().data() + s_offset, m_bsplines.nbasis());
    matrix->solve_inplace(bcoef_section);

    if constexpr (BcXmin == BoundCond::PERIODIC && s_offset != 0) {
        for (int i(0); i < s_offset; ++i) {
            spline(i) = spline(m_bsplines.nbasis() + i);
        }
        for (int i(s_offset); i < m_bsplines.degree(); ++i) {
            spline(m_bsplines.nbasis() + i) = spline(i);
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
    int const n_interp_pts = m_bsplines.nbasis() - s_nbc_xmin - s_nbc_xmax;

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        double const shift(!s_odd ? 0.5 : 0.0);
        UniformMesh<tag_type> mesh(m_bsplines.rmin() + shift * m_dx, m_dx);
        m_interpolation_domain
                = std::make_unique<interpolation_domain_type>(mesh, MCoord<tag_type>(n_interp_pts));
    } else {
        std::vector<double> interp_pts(n_interp_pts);

        int n_iknots = n_interp_pts + m_bsplines.degree() - 1;
        std::vector<int> iknots(n_iknots);
        int i(0);

        // Additional knots near x=xmin
        int n_to_fill_min(m_bsplines.degree() - s_nbc_xmin - 1);
        for (; i < n_to_fill_min; ++i) {
            if constexpr (BcXmin == BoundCond::GREVILLE)
                iknots[i] = 0;
            if constexpr (BcXmin == BoundCond::HERMITE)
                iknots[i] = -n_to_fill_min + i;
        }

        // Knots inside the domain
        for (int j(0); j < m_bsplines.ncells() + 1; ++i, ++j) {
            iknots[i] = j;
        }

        // Additional knots near x=xmax
        for (int j(1); i < n_iknots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                iknots[i] = m_bsplines.ncells();
            if constexpr (BcXmax == BoundCond::HERMITE)
                iknots[i] = m_bsplines.ncells() + j;
        }

        for (int j(0); j < n_interp_pts; ++j) {
            int isum(sum(iknots.data() + j, m_bsplines.degree()));
            interp_pts[j] = m_bsplines.rmin() + m_dx * isum / m_bsplines.degree();
        }

        // Non-periodic case, odd degree: fix round-off issues
        if constexpr (s_odd) {
            interp_pts[0] = m_bsplines.rmin();
            interp_pts[n_interp_pts - 1] = m_bsplines.rmax();
        }
        NonUniformMesh<tag_type> mesh(interp_pts, MCoord<tag_type>(0));
        m_interpolation_domain = std::make_unique<
                interpolation_domain_type>(mesh, MCoord<tag_type>(interp_pts.size()));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, BcXmin, BcXmax>::compute_interpolation_points_non_uniform()
{
    int n_interp_pts = m_bsplines.nbasis() - s_nbc_xmin - s_nbc_xmax;
    std::vector<double> interp_pts(n_interp_pts);

    int n_temp_knots(n_interp_pts - 1 + m_bsplines.degree());
    double temp_knots[n_temp_knots];

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        for (int i(0); i < n_interp_pts - 1 + m_bsplines.degree(); ++i) {
            temp_knots[i] = m_bsplines.get_knot(1 - m_bsplines.degree() + s_offset + i);
        }
    } else {
        int i(0);
        int n_start_pts(m_bsplines.degree() - s_nbc_xmin - 1);

        // Initialise knots relevant to the xmin boundary condition
        for (; i < n_start_pts; ++i) {
            // As xmin_bc is a const variable the compiler should optimize
            // for(if..else..) to if(for..)else(for...)
            if constexpr (BcXmin == BoundCond::GREVILLE)
                temp_knots[i] = m_bsplines.get_knot(0);
            if constexpr (BcXmin == BoundCond::HERMITE)
                temp_knots[i] = 2.0 * m_bsplines.get_knot(0) - m_bsplines.get_knot(n_start_pts - i);
        }

        // Initialise central knots
        for (int j(0); j < m_bsplines.npoints(); ++i, ++j) {
            temp_knots[i] = m_bsplines.get_knot(j);
        }

        // Initialise knots relevant to the xmax boundary condition
        for (int j(0); i < n_temp_knots; ++i, ++j) {
            if constexpr (BcXmax == BoundCond::GREVILLE)
                temp_knots[i] = m_bsplines.get_knot(m_bsplines.ncells());
            if constexpr (BcXmax == BoundCond::HERMITE)
                temp_knots[i] = 2.0 * m_bsplines.get_knot(m_bsplines.ncells())
                                - m_bsplines.get_knot(m_bsplines.ncells() - 1 - j);
        }
    }

    // Compute interpolation points using Greville-style averaging
    double inv_deg = 1.0 / m_bsplines.degree();
    for (int i(0); i < n_interp_pts; ++i) {
        interp_pts[i] = sum(temp_knots + i, m_bsplines.degree()) * inv_deg;
    }

    // Periodic case: apply periodic BCs to interpolation points
    if constexpr (BcXmin == BoundCond::PERIODIC) {
        double zone_width(m_bsplines.rmax() - m_bsplines.rmin());
        for (int i(0); i < n_interp_pts; ++i) {
            interp_pts[i] = modulo(interp_pts[i] - s_nbc_xmin, zone_width) + m_bsplines.rmin();
        }
    }
    // Non-periodic case, odd degree: fix round-off issues
    else {
        if constexpr (s_odd) {
            interp_pts[0] = m_bsplines.rmin();
            interp_pts[n_interp_pts - 1] = m_bsplines.rmax();
        }
    }

    NonUniformMesh<tag_type> mesh(interp_pts, MCoord<tag_type>(0));
    m_interpolation_domain = std::make_unique<
            interpolation_domain_type>(mesh, MCoord<tag_type>(interp_pts.size()));
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
        upper_block_size = (m_bsplines.degree()) / 2;
        break;
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = m_bsplines.degree() - 1;
        break;
    default:
        break; // TODO: throw error
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = (m_bsplines.degree()) / 2;
        break;
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = m_bsplines.degree() - 1;
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
        upper_block_size = (m_bsplines.degree() + 1) / 2;
        break;
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin + 1;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = m_bsplines.degree() - 1;
        break;
    default:
        break; // TODO: throw error
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = (m_bsplines.degree() + 1) / 2;
        break;
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax + 1;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = m_bsplines.degree() - 1;
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
    if (m_bsplines.is_uniform()) {
        upper_band_width = m_bsplines.degree() / 2;
    } else {
        upper_band_width = (m_bsplines.degree() + 1) / 2;
    }

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        matrix = Matrix::make_new_periodic_banded(
                m_bsplines.nbasis(),
                upper_band_width,
                upper_band_width,
                m_bsplines.is_uniform());
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                m_bsplines.nbasis(),
                upper_band_width,
                upper_band_width,
                m_bsplines.is_uniform(),
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
        double derivs_ptr[(m_bsplines.degree() / 2 + 1) * (m_bsplines.degree() + 1)];
        DSpan2D derivs(derivs_ptr, m_bsplines.degree() + 1, m_bsplines.degree() / 2 + 1);
        m_bsplines.eval_basis_and_n_derivs(m_bsplines.rmin(), s_nbc_xmin, derivs, jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i(0); i < m_bsplines.degree() + 1; ++i) {
            for (int j(1); j < m_bsplines.degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        // iterate only to deg as last bspline is 0
        for (int j(0); j < s_nbc_xmin; ++j) {
            for (int i(0); i < m_bsplines.degree(); ++i) {
                // Elements are set in Fortran order as they are LAPACK input
                matrix->set_element(i, j, derivs(i, s_nbc_xmin - j - 1 + s_odd));
            }
        }
    }

    // Interpolation points
    double values_ptr[m_bsplines.degree() + 1];
    DSpan1D values(values_ptr, m_bsplines.degree() + 1);
    for (int i(0); i < m_bsplines.nbasis() - s_nbc_xmin - s_nbc_xmax; ++i) {
        m_bsplines.eval_basis(m_interpolation_domain->to_real(i), values, jmin);
        for (int s(0); s < m_bsplines.degree() + 1; ++s) {
            int j = modulo(jmin - s_offset + s, (int)m_bsplines.nbasis());
            matrix->set_element(j, i + s_nbc_xmin, values(s));
        }
    }

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == BoundCond::HERMITE) {
        double derivs_ptr[(m_bsplines.degree() / 2 + 1) * (m_bsplines.degree() + 1)];
        DSpan2D derivs(derivs_ptr, m_bsplines.degree() + 1, m_bsplines.degree() / 2 + 1);

        m_bsplines.eval_basis_and_n_derivs(m_bsplines.rmax(), s_nbc_xmax, derivs, jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i(0); i < m_bsplines.degree() + 1; ++i) {
            for (int j(1); j < m_bsplines.degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
            }
        }

        int i0(m_bsplines.nbasis() - m_bsplines.degree());
        int j0(m_bsplines.nbasis() - s_nbc_xmax);
        for (int i(0); i < m_bsplines.degree(); ++i) {
            for (int j(0); j < s_nbc_xmax; ++j) {
                matrix->set_element(i0 + i, j0 + j, derivs(i + 1, j + s_odd));
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                                 Static functions *
 ************************************************************************************/

template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
int SplineBuilder<BSplines, BcXmin, BcXmax>::compute_num_cells(
        int degree,
        BoundCond xmin_bc,
        BoundCond xmax_bc,
        int nipts)
{
    assert(degree > 0);
    // TODO: xmin in allowed_bcs
    // TODO: xmax in allowed_bcs

    if constexpr ((BcXmin == BoundCond::PERIODIC) != (BcXmax == BoundCond::PERIODIC)) {
        std::cerr << "Incompatible BCs" << std::endl;
        // TODO: raise error
        return -1;
    }

    if constexpr (BcXmin == BoundCond::PERIODIC) {
        return nipts;
    } else {
        int nbc_xmin, nbc_xmax;
        if constexpr (BcXmin == BoundCond::HERMITE)
            nbc_xmin = degree / 2;
        else
            nbc_xmin = 0;

        if constexpr (BcXmax == BoundCond::HERMITE)
            nbc_xmax = degree / 2;
        else
            nbc_xmax = 0;

        return nipts + nbc_xmin + nbc_xmax - degree;
    }
}
