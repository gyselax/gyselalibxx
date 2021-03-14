#include <cassert>
#include <iostream>

#include <matrix_banded.h>
#include <matrix_dense.h>
#include <matrix_periodic_banded.h>

#include "spline_interpolator_1d.h"

std::array<BoundaryCondition, 3> Spline_interpolator_1D::allowed_bcs
        = {BoundaryCondition::sll_p_periodic,
           BoundaryCondition::sll_p_hermite,
           BoundaryCondition::sll_p_greville};

Spline_interpolator_1D::Spline_interpolator_1D(
        const BSplines& bspl,
        BoundaryCondition xmin_bc,
        BoundaryCondition xmax_bc)
    : xmin_bc(xmin_bc)
    , xmax_bc(xmax_bc)
    , nbc_xmin(xmin_bc == BoundaryCondition::sll_p_hermite ? bspl.degree / 2 : 0)
    , nbc_xmax(xmax_bc == BoundaryCondition::sll_p_hermite ? bspl.degree / 2 : 0)
    , bspl(bspl)
    , odd(bspl.degree % 2)
    , offset(bspl.periodic ? bspl.degree / 2 : 0)
    , dx((bspl.xmax - bspl.xmin) / bspl.ncells)
    , matrix(nullptr)
{
    int lower_block_size, upper_block_size;
    constructor_sanity_checks();
    if (bspl.uniform) {
        compute_interpolation_points_uniform();
        compute_block_sizes_uniform(lower_block_size, upper_block_size);
    } else {
        compute_interpolation_points_non_uniform();
        compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
    }
    allocate_matrix(lower_block_size, upper_block_size);
}

//-------------------------------------------------------------------------------------------------

inline void Spline_interpolator_1D::constructor_sanity_checks()
{
    assert(xmin_bc == allowed_bcs[0] || xmin_bc == allowed_bcs[1] || xmin_bc == allowed_bcs[2]);
    assert(xmax_bc == allowed_bcs[0] || xmax_bc == allowed_bcs[1] || xmax_bc == allowed_bcs[2]);
    if (bspl.periodic) {
        assert(xmin_bc == BoundaryCondition::sll_p_periodic);
        assert(xmax_bc == BoundaryCondition::sll_p_periodic);
    }

    assert(not bspl.radial);
}

//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------

const mdspan_1d& Spline_interpolator_1D::get_interp_points() const
{
    return interp_pts;
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

void Spline_interpolator_1D::compute_interpolant_degree1(Spline_1D& spline, const mdspan_1d& vals)
        const
{
    for (int i(0); i < bspl.nbasis; ++i) {
        spline.bcoef(i) = vals(i);
    }
    if (bspl.periodic) {
        spline.bcoef(bspl.nbasis) = spline.bcoef(0);
    }
    return;
}

//-------------------------------------------------------------------------------------------------

void Spline_interpolator_1D::compute_interpolant(
        Spline_1D& spline,
        const mdspan_1d& vals,
        const mdspan_1d* derivs_xmin,
        const mdspan_1d* derivs_xmax) const
{
    assert(vals.extent(0) == bspl.nbasis - nbc_xmin - nbc_xmax);
    assert(spline.belongs_to_space(bspl));
    // TODO: LOG Errors
    if (bspl.degree == 1)
        return compute_interpolant_degree1(spline, vals);

    assert((xmin_bc == BoundaryCondition::sll_p_hermite)
           != (derivs_xmin == nullptr || derivs_xmin->extent(0) == 0));
    assert((xmax_bc == BoundaryCondition::sll_p_hermite)
           != (derivs_xmax == nullptr || derivs_xmax->extent(0) == 0));

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if (xmin_bc == BoundaryCondition::sll_p_hermite) {
        for (int i(nbc_xmin); i > 0; --i) {
            spline.bcoef(nbc_xmin - i) = (*derivs_xmin)(i - 1) * ipow(dx, i + odd - 1);
        }
    }
    for (int i(nbc_xmin); i < nbc_xmin + offset; ++i) {
        spline.bcoef(i) = 0.0;
    }

    for (int i(0); i < interp_pts.extent(0); ++i) {
        spline.bcoef(nbc_xmin + i + offset) = vals(i);
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if (xmax_bc == BoundaryCondition::sll_p_hermite) {
        for (int i(0); i < nbc_xmax; ++i) {
            spline.bcoef(bspl.nbasis - nbc_xmax + i) = (*derivs_xmax)(i)*ipow(dx, i + odd);
        }
    }

    mdspan_1d bcoef_section(spline.bcoef_ptr.get() + offset, bspl.nbasis);
    matrix->solve_inplace(bcoef_section);

    if (xmin_bc == BoundaryCondition::sll_p_periodic and offset != 0) {
        for (int i(0); i < offset; ++i) {
            spline.bcoef(i) = spline.bcoef(bspl.nbasis + i);
        }
        for (int i(offset); i < bspl.degree; ++i) {
            spline.bcoef(bspl.nbasis + i) = spline.bcoef(i);
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                  Compute interpolation points functions *
 ************************************************************************************/

void Spline_interpolator_1D::compute_interpolation_points_uniform()
{
    int n_interp_pts = bspl.nbasis - nbc_xmin - nbc_xmax;
    interp_pts_ptr = std::make_unique<double[]>(n_interp_pts);
    interp_pts = mdspan_1d(interp_pts_ptr.get(), n_interp_pts);

    if (xmin_bc == BoundaryCondition::sll_p_periodic) {
        double shift(odd == 0 ? 0.5 : 0.0);
        for (int i(0); i < n_interp_pts; ++i) {
            interp_pts(i) = bspl.xmin + (i + shift) * dx;
        }
    } else {
        int n_iknots = n_interp_pts + bspl.degree - 1;
        int iknots[n_iknots];
        int i(0);

        // Additional knots near x=xmin
        int n_to_fill_min(bspl.degree - nbc_xmin - 1);
        for (; i < n_to_fill_min; ++i) {
            if (xmin_bc == BoundaryCondition::sll_p_greville)
                iknots[i] = 0;
            if (xmin_bc == BoundaryCondition::sll_p_hermite)
                iknots[i] = -n_to_fill_min + i;
        }

        // Knots inside the domain
        for (int j(0); j < bspl.ncells + 1; ++i, ++j) {
            iknots[i] = j;
        }

        // Additional knots near x=xmax
        for (int j(1); i < n_iknots; ++i, ++j) {
            if (xmax_bc == BoundaryCondition::sll_p_greville)
                iknots[i] = bspl.ncells;
            if (xmax_bc == BoundaryCondition::sll_p_hermite)
                iknots[i] = bspl.ncells + j;
        }

        for (int j(0); j < n_interp_pts; ++j) {
            int isum(sum(iknots + j, bspl.degree));
            interp_pts(j) = bspl.xmin + dx * isum / bspl.degree;
        }

        // Non-periodic case, odd degree: fix round-off issues
        if (odd == 1) {
            interp_pts(0) = bspl.xmin;
            interp_pts(n_interp_pts - 1) = bspl.xmax;
        }
    }
}

//-------------------------------------------------------------------------------------------------

void Spline_interpolator_1D::compute_interpolation_points_non_uniform()
{
    const BSplines_non_uniform& bspl_nu = static_cast<const BSplines_non_uniform&>(bspl);
    int n_interp_pts = bspl_nu.nbasis - nbc_xmin - nbc_xmax;
    interp_pts_ptr = std::make_unique<double[]>(n_interp_pts);
    interp_pts = mdspan_1d(interp_pts_ptr.get(), n_interp_pts);

    int n_temp_knots(n_interp_pts - 1 + bspl_nu.degree);
    double temp_knots[n_temp_knots];

    if (xmin_bc == BoundaryCondition::sll_p_periodic) {
        for (int i(0); i < n_interp_pts - 1 + bspl_nu.degree; ++i) {
            temp_knots[i] = bspl_nu.get_knot(1 - bspl_nu.degree + offset + i);
        }
    } else {
        int i(0);
        int n_start_pts(bspl_nu.degree - nbc_xmin - 1);

        // Initialise knots relevant to the xmin boundary condition
        for (; i < n_start_pts; ++i) {
            // As xmin_bc is a const variable the compiler should optimize
            // for(if..else..) to if(for..)else(for...)
            if (xmin_bc == BoundaryCondition::sll_p_greville)
                temp_knots[i] = bspl_nu.get_knot(0);
            if (xmin_bc == BoundaryCondition::sll_p_hermite)
                temp_knots[i] = 2.0 * bspl_nu.get_knot(0) - bspl_nu.get_knot(n_start_pts - i);
        }

        // Initialise central knots
        for (int j(0); j < bspl_nu.npoints; ++i, ++j) {
            temp_knots[i] = bspl_nu.get_knot(j);
        }

        // Initialise knots relevant to the xmax boundary condition
        for (int j(0); i < n_temp_knots; ++i, ++j) {
            if (xmax_bc == BoundaryCondition::sll_p_greville)
                temp_knots[i] = bspl_nu.get_knot(bspl_nu.ncells);
            if (xmax_bc == BoundaryCondition::sll_p_hermite)
                temp_knots[i] = 2.0 * bspl_nu.get_knot(bspl_nu.ncells)
                              - bspl_nu.get_knot(bspl_nu.ncells - 1 - j);
        }
    }

    // Compute interpolation points using Greville-style averaging
    double inv_deg = 1.0 / bspl.degree;
    for (int i(0); i < n_interp_pts; ++i) {
        interp_pts(i) = sum(temp_knots + i, bspl.degree) * inv_deg;
    }

    // Periodic case: apply periodic BCs to interpolation points
    if (xmin_bc == BoundaryCondition::sll_p_periodic) {
        double zone_width(bspl.xmax - bspl.xmin);
        for (int i(0); i < n_interp_pts; ++i) {
            interp_pts(i) = modulo(interp_pts(i) - nbc_xmin, zone_width) + bspl.xmin;
        }
    }
    // Non-periodic case, odd degree: fix round-off issues
    else if (odd == 1) {
        interp_pts(0) = bspl.xmin;
        interp_pts(n_interp_pts - 1) = bspl.xmax;
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Compute num diags functions *
 ************************************************************************************/

void Spline_interpolator_1D::compute_block_sizes_uniform(
        int& lower_block_size,
        int& upper_block_size) const
{
    switch (xmin_bc) {
    case BoundaryCondition::sll_p_periodic:
        upper_block_size = (bspl.degree) / 2;
        break;
    case BoundaryCondition::sll_p_hermite:
        upper_block_size = nbc_xmin;
        break;
    case BoundaryCondition::sll_p_greville:
        upper_block_size = bspl.degree - 1;
        break;
    default:
        break; // TODO: throw error
    }
    switch (xmax_bc) {
    case BoundaryCondition::sll_p_periodic:
        lower_block_size = (bspl.degree) / 2;
        break;
    case BoundaryCondition::sll_p_hermite:
        lower_block_size = nbc_xmax;
        break;
    case BoundaryCondition::sll_p_greville:
        lower_block_size = bspl.degree - 1;
        break;
    default:
        break; // TODO: throw error
    }
}

//-------------------------------------------------------------------------------------------------

void Spline_interpolator_1D::compute_block_sizes_non_uniform(
        int& lower_block_size,
        int& upper_block_size) const
{
    switch (xmin_bc) {
    case BoundaryCondition::sll_p_periodic:
        upper_block_size = (bspl.degree + 1) / 2;
        break;
    case BoundaryCondition::sll_p_hermite:
        upper_block_size = nbc_xmin + 1;
        break;
    case BoundaryCondition::sll_p_greville:
        upper_block_size = bspl.degree - 1;
        break;
    default:
        break; // TODO: throw error
    }
    switch (xmax_bc) {
    case BoundaryCondition::sll_p_periodic:
        lower_block_size = (bspl.degree + 1) / 2;
        break;
    case BoundaryCondition::sll_p_hermite:
        lower_block_size = nbc_xmax + 1;
        break;
    case BoundaryCondition::sll_p_greville:
        lower_block_size = bspl.degree - 1;
        break;
    default:
        break; // TODO: throw error
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Initailise matrix functions *
 ************************************************************************************/

void Spline_interpolator_1D::allocate_matrix(int lower_block_size, int upper_block_size)
{
    // Special case: linear spline
    // No need for matrix assembly
    if (bspl.degree == 1)
        return;

    int upper_band_width;
    if (bspl.uniform) {
        upper_band_width = bspl.degree / 2;
    } else {
        upper_band_width = (bspl.degree + 1) / 2;
    }

    if (xmin_bc == BoundaryCondition::sll_p_periodic) {
        matrix = Matrix::make_new_periodic_banded(
                bspl.nbasis,
                upper_band_width,
                upper_band_width,
                bspl.uniform);
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                bspl.nbasis,
                upper_band_width,
                upper_band_width,
                bspl.uniform,
                upper_block_size,
                lower_block_size);
    }

    build_matrix_system();

    matrix->factorize();
}

//-------------------------------------------------------------------------------------------------

void Spline_interpolator_1D::build_matrix_system()
{
    int jmin;

    // Hermite boundary conditions at xmin, if any
    if (xmin_bc == BoundaryCondition::sll_p_hermite) {
        double derivs_ptr[(bspl.degree / 2 + 1) * (bspl.degree + 1)];
        mdspan_2d derivs(derivs_ptr, bspl.degree + 1, bspl.degree / 2 + 1);
        bspl.eval_basis_and_n_derivs(bspl.xmin, nbc_xmin, derivs, jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i(0); i < bspl.degree + 1; ++i) {
            for (int j(1); j < bspl.degree / 2 + 1; ++j) {
                derivs(i, j) *= ipow(dx, j);
            }
        }

        // iterate only to deg as last bspline is 0
        for (int j(0); j < nbc_xmin; ++j) {
            for (int i(0); i < bspl.degree; ++i) {
                // Elements are set in Fortran order as they are LAPACK input
                matrix->set_element(i, j, derivs(i, nbc_xmin - j - 1 + odd));
            }
        }
    }

    // Interpolation points
    double values_ptr[bspl.degree + 1];
    mdspan_1d values(values_ptr, bspl.degree + 1);
    for (int i(0); i < bspl.nbasis - nbc_xmin - nbc_xmax; ++i) {
        bspl.eval_basis(interp_pts(i), values, jmin);
        for (int s(0); s < bspl.degree + 1; ++s) {
            int j = modulo(jmin - offset + s, bspl.nbasis);
            matrix->set_element(j, i + nbc_xmin, values(s));
        }
    }

    // Hermite boundary conditions at xmax, if any
    if (xmax_bc == BoundaryCondition::sll_p_hermite) {
        double derivs_ptr[(bspl.degree / 2 + 1) * (bspl.degree + 1)];
        mdspan_2d derivs(derivs_ptr, bspl.degree + 1, bspl.degree / 2 + 1);

        bspl.eval_basis_and_n_derivs(bspl.xmax, nbc_xmax, derivs, jmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (int i(0); i < bspl.degree + 1; ++i) {
            for (int j(1); j < bspl.degree / 2 + 1; ++j) {
                derivs(i, j) *= ipow(dx, j);
            }
        }

        int i0(bspl.nbasis - bspl.degree);
        int j0(bspl.nbasis - nbc_xmax);
        for (int i(0); i < bspl.degree; ++i) {
            for (int j(0); j < nbc_xmax; ++j) {
                matrix->set_element(i0 + i, j0 + j, derivs(i + 1, j + odd));
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                                 Static functions *
 ************************************************************************************/

int Spline_interpolator_1D::compute_num_cells(
        int degree,
        BoundaryCondition xmin_bc,
        BoundaryCondition xmax_bc,
        int nipts)
{
    assert(degree > 0);
    // TODO: xmin in allowed_bcs
    // TODO: xmax in allowed_bcs

    if ((xmin_bc == BoundaryCondition::sll_p_periodic)
        != (xmax_bc == BoundaryCondition::sll_p_periodic)) {
        std::cerr << "Incompatible BCs" << std::endl;
        // TODO: raise error
        return -1;
    }

    if (xmin_bc == BoundaryCondition::sll_p_periodic) {
        return nipts;
    } else {
        int nbc_xmin, nbc_xmax;
        if (xmin_bc == BoundaryCondition::sll_p_hermite)
            nbc_xmin = degree / 2;
        else
            nbc_xmin = 0;

        if (xmax_bc == BoundaryCondition::sll_p_hermite)
            nbc_xmax = degree / 2;
        else
            nbc_xmax = 0;

        return nipts + nbc_xmin + nbc_xmax - degree;
    }
}
