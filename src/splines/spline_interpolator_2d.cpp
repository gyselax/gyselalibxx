#include <cassert>

#include "spline_interpolator_2d.h"

/******************************************************************************
 * @brief      Initialize a 2D tensor product spline interpolator object
 * @param[out] self     2D tensor product spline interpolator
 * @param[in]  bspl1    B-splines (basis) along x1 direction
 * @param[in]  bspl2    B-splines (basis) along x2 direction
 * @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
 * @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
 ******************************************************************************/
Spline_interpolator_2D::Spline_interpolator_2D(
        std::array<std::unique_ptr<const BSplines>, 2> bspl,
        std::array<BoundCond, 2> xmin_bc,
        std::array<BoundCond, 2> xmax_bc)
    : bspl(std::move(bspl))
    , interp_1d(
              {Spline_interpolator_1D(*bspl[0], xmin_bc[0], xmax_bc[0]),
               Spline_interpolator_1D(*bspl[1], xmin_bc[1], xmax_bc[1])})
    , spline_1d({Spline_1D(*bspl[0]), Spline_1D(*bspl[1])})
    , xmin_bc(xmin_bc)
    , xmax_bc(xmax_bc)
    , nbc_xmin({interp_1d[0].nbc_xmin, interp_1d[1].nbc_xmin})
    , nbc_xmax({interp_1d[0].nbc_xmax, interp_1d[1].nbc_xmax})
{
}

/******************************************************************************
 * @brief      Get coordinates of interpolation points (2D tensor grid)
 * @param[in]  self  2D tensor product spline interpolator
 * @param[out] tau1  x1 coordinates of interpolation points
 * @param[out] tau2  x2 coordinates of interpolation points
 ******************************************************************************/
std::array<const DSpan1D, 2> Spline_interpolator_2D::get_interp_points() const
{
    return {interp_1d[0].get_interp_points(), interp_1d[1].get_interp_points()};
}

/******************************************************************************
 * @brief      Calculate number of cells from number of interpolation points
 * @details    Important for parallelization: for given spline degree and BCs,
 *             calculate the numbers of grid cells along x1 and x2 that yield
 *             the desired number of interpolation points along x1 and x2
 *
 * @param[in]  degree   spline degrees along x1 and x2
 * @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
 * @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
 * @param[in]  nipts    desired number of interpolation points along x1 and x2
 * @param[out] ncells   calculated number of grid cells along x1 and x2
 ******************************************************************************/
std::array<int, 2> Spline_interpolator_2D::compute_num_cells(
        std::array<int, 2> degree,
        std::array<BoundCond, 2> xmin,
        std::array<BoundCond, 2> xmax,
        std::array<int, 2> nipts)
{
    int n0(Spline_interpolator_1D::compute_num_cells(degree[0], xmin[0], xmax[0], nipts[0]));
    int n1(Spline_interpolator_1D::compute_num_cells(degree[1], xmin[1], xmax[1], nipts[1]));
    return {n0, n1};
}

/******************************************************************************
 * @brief        Compute interpolating 2D spline
 * @details      Compute coefficients of 2D tensor product spline that
 *               interpolates function values on grid. If Hermite BCs are used,
 *               function derivatives at appropriate boundaries are also
 *needed.
 *
 * @param[inout] self           2D tensor product spline interpolator
 * @param[inout] spline         2D tensor product spline
 * @param[in]    gtau           function values of interpolation points
 * @param[in]    boundary_data  (optional) structure with boundary conditions
 ******************************************************************************/
void Spline_interpolator_2D::compute_interpolant(
        Spline_2D const& spline,
        DSpan2D const& vals,
        Boundary_data_2d boundary_data) const
{
    int dim_0_size = (bspl[0]->nbasis - nbc_xmin[0] - nbc_xmax[0]);
    int dim_1_size = (bspl[1]->nbasis - nbc_xmin[1] - nbc_xmax[1]);
    // TODO fix assert
    if (nbc_xmin[0] > 0) {
        assert(boundary_data.derivs_x1_min != nullptr);
        assert(boundary_data.derivs_x1_min->extent(0) == dim_1_size);
        assert(boundary_data.derivs_x1_min->extent(1) == nbc_xmin[0]);
        for (int i(0); i < nbc_xmin[0]; ++i) {
            for (int j(0); j < dim_1_size; ++j) {
                spline.bcoef(i, j + nbc_xmin[1]) = (*boundary_data.derivs_x1_min)(j, i);
            }
        }
    }
    if (nbc_xmax[0] > 0) {
        assert(boundary_data.derivs_x1_max != nullptr);
        assert(boundary_data.derivs_x1_max->extent(0) == dim_1_size);
        assert(boundary_data.derivs_x1_max->extent(1) == nbc_xmax[0]);
        for (int i(0); i < nbc_xmax[0]; ++i) {
            for (int j(0); j < dim_1_size; ++j) {
                spline.bcoef(i + bspl[0]->nbasis - nbc_xmax[0], j + nbc_xmin[1])
                        = (*boundary_data.derivs_x1_max)(j, i);
            }
        }
    }
    if (nbc_xmin[1] > 0) {
        assert(boundary_data.derivs_x2_min != nullptr);
        assert(boundary_data.derivs_x2_min->extent(0) == dim_0_size);
        assert(boundary_data.derivs_x2_min->extent(1) == nbc_xmin[1]);
        for (int i(0); i < dim_0_size; ++i) {
            for (int j(0); j < nbc_xmin[1]; ++j) {
                spline.bcoef(i + nbc_xmin[0], j) = (*boundary_data.derivs_x2_min)(i, j);
            }
        }
    }
    if (nbc_xmax[1] > 0) {
        assert(boundary_data.derivs_x2_max != nullptr);
        assert(boundary_data.derivs_x2_max->extent(0) == dim_0_size);
        assert(boundary_data.derivs_x2_max->extent(1) == nbc_xmax[1]);
        for (int i(0); i < dim_0_size; ++i) {
            for (int j(0); j < nbc_xmax[1]; ++j) {
                spline.bcoef(i + nbc_xmin[0], j + bspl[1]->nbasis - nbc_xmax[1])
                        = (*boundary_data.derivs_x2_max)(i, j);
            }
        }
    }
    if (nbc_xmin[0] > 0 and nbc_xmin[1] > 0) {
        assert(boundary_data.mixed_derivs_a != nullptr);
        assert(boundary_data.mixed_derivs_a->extent(0) == nbc_xmin[0]);
        assert(boundary_data.mixed_derivs_a->extent(1) == nbc_xmin[1]);
        for (int i(0); i < nbc_xmin[0]; ++i) {
            for (int j(0); j < nbc_xmin[1]; ++j) {
                spline.bcoef(i, j) = (*boundary_data.mixed_derivs_a)(i, j);
            }
        }
    }
    if (nbc_xmax[0] > 0 and nbc_xmin[1] > 0) {
        assert(boundary_data.mixed_derivs_b != nullptr);
        assert(boundary_data.mixed_derivs_b->extent(0) == nbc_xmax[0]);
        assert(boundary_data.mixed_derivs_b->extent(1) == nbc_xmin[1]);
        for (int i(0); i < nbc_xmax[0]; ++i) {
            for (int j(0); j < nbc_xmin[1]; ++j) {
                spline.bcoef(i + bspl[0]->nbasis - nbc_xmax[0], j)
                        = (*boundary_data.mixed_derivs_b)(i, j);
            }
        }
    }
    if (nbc_xmin[0] > 0 and nbc_xmax[1] > 0) {
        assert(boundary_data.mixed_derivs_c != nullptr);
        assert(boundary_data.mixed_derivs_c->extent(0) == nbc_xmin[0]);
        assert(boundary_data.mixed_derivs_c->extent(1) == nbc_xmax[1]);
        for (int i(0); i < nbc_xmin[0]; ++i) {
            for (int j(0); j < nbc_xmax[1]; ++j) {
                spline.bcoef(i, j + bspl[1]->nbasis - nbc_xmax[1])
                        = (*boundary_data.mixed_derivs_c)(i, j);
            }
        }
    }
    if (nbc_xmax[0] > 0 and nbc_xmax[1] > 0) {
        assert(boundary_data.mixed_derivs_d != nullptr);
        assert(boundary_data.mixed_derivs_d->extent(0) == nbc_xmax[0]);
        assert(boundary_data.mixed_derivs_d->extent(1) == nbc_xmax[1]);
        for (int i(0); i < nbc_xmax[0]; ++i) {
            for (int j(0); j < nbc_xmax[1]; ++j) {
                spline.bcoef(i + bspl[0]->nbasis - nbc_xmax[0], j + bspl[1]->nbasis - nbc_xmax[1])
                        = (*boundary_data.mixed_derivs_d)(i, j);
            }
        }
    }

    compute_interpolant_boundary_done(spline, vals);
}

/******************************************************************************
 * @brief        Compute interpolating 2D spline
 * @details      Compute coefficients of 2D tensor product spline that
 *               interpolates function values on grid. If Hermite BCs are used,
 *               function derivatives at appropriate boundaries are also
 *needed.
 *
 * @param[inout] self           2D tensor product spline interpolator
 * @param[inout] spline         2D tensor product spline
 * @param[in]    gtau           function values of interpolation points
 * @param[in]    boundary_data  (optional) structure with boundary conditions
 ******************************************************************************/
void Spline_interpolator_2D::compute_interpolant(Spline_2D const& spline, DSpan2D const& vals) const
{
    assert(xmin_bc[0] != BoundCond::HERMITE);
    assert(xmax_bc[0] != BoundCond::HERMITE);
    assert(xmin_bc[1] != BoundCond::HERMITE);
    assert(xmax_bc[1] != BoundCond::HERMITE);
    compute_interpolant_boundary_done(spline, vals);
}

/******************************************************************************
 * @brief        Compute interpolating 2D spline
 * @details      Compute coefficients of 2D tensor product spline that
 *               interpolates function values on grid. If Hermite BCs are used,
 *               function derivatives at appropriate boundaries are also
 *needed.
 *
 * @param[inout] self           2D tensor product spline interpolator
 * @param[inout] spline         2D tensor product spline
 * @param[in]    gtau           function values of interpolation points
 * @param[in]    boundary_data  (optional) structure with boundary conditions
 ******************************************************************************/
void Spline_interpolator_2D::compute_interpolant_boundary_done(
        Spline_2D const& spline,
        DSpan2D const& vals) const
{
    assert(vals.extent(0) == (bspl[0]->nbasis - nbc_xmin[0] - nbc_xmax[0]));
    assert(vals.extent(1) == (bspl[1]->nbasis - nbc_xmin[1] - nbc_xmax[1]));
    assert(spline.belongs_to_space(*bspl[0], *bspl[1]));

    int dim_0_size = (bspl[0]->nbasis - nbc_xmin[0] - nbc_xmax[0]);
    int dim_1_size = (bspl[1]->nbasis - nbc_xmin[1] - nbc_xmax[1]);

    // Copy interpolation data onto w array
    for (int i(0); i < vals.extent(0); ++i) {
        for (int j(0); j < vals.extent(1); ++j) {
            spline.bcoef(nbc_xmin[0] + i, nbc_xmin[1] + j) = vals(i, j);
        }
    }

    std::array<Spline_1D, 2> spline_1d({Spline_1D(*bspl[0]), Spline_1D(*bspl[1])});
    double t_storage_ptr[spline.bcoef.extent(0) * spline.bcoef.extent(1)];
    DSpan2D t_storage(t_storage_ptr, spline.bcoef.extent(1), spline.bcoef.extent(0));
    // Cycle over x1 position (or order of x1-derivative at boundary)
    // and interpolate f along x2 direction. Store coefficients in bcoef
    {
        int i(0);
        for (; i < bspl[0]->nbasis; ++i) {
            DSpan1D values(&spline.bcoef(i, nbc_xmin[1]), dim_1_size);
            DSpan1D derivs_xmin(&spline.bcoef(i, 0), nbc_xmin[1]);
            DSpan1D derivs_xmax(&spline.bcoef(i, bspl[1]->nbasis - nbc_xmax[1]), nbc_xmax[1]);
            interp_1d[1].compute_interpolant(spline_1d[1], values, &derivs_xmin, &derivs_xmax);

            int j(0);
            for (; j < bspl[1]->nbasis; ++j) {
                t_storage(j, i) = spline_1d[1].bcoef(j);
            }
            for (; j < spline.bcoef.extent(1); ++j) {
                t_storage(j, i) = spline.bcoef(i, j);
            }
        }
        for (; i < spline.bcoef.extent(0); ++i) {
            for (int j(0); j < spline.bcoef.extent(1); ++j) {
                t_storage(j, i) = spline.bcoef(i, j);
            }
        }
    }

    // Cycle over x2 position (or order of x2-derivative at boundary)
    // and interpolate f along x1 direction. Store coefficients in bwork
    {
        int j(0);
        for (; j < bspl[1]->nbasis; ++j) {
            DSpan1D values(&t_storage(j, nbc_xmin[0]), dim_0_size);
            DSpan1D derivs_xmin(&t_storage(j, 0), nbc_xmin[0]);
            DSpan1D derivs_xmax(&t_storage(j, bspl[0]->nbasis - nbc_xmax[0]), nbc_xmax[0]);
            interp_1d[0].compute_interpolant(spline_1d[0], values, &derivs_xmin, &derivs_xmax);

            int i(0);
            for (; i < bspl[0]->nbasis; ++i) {
                spline.bcoef(i, j) = spline_1d[0].bcoef(i);
            }
            for (; i < spline.bcoef.extent(0); ++i) {
                spline.bcoef(i, j) = t_storage(j, i);
            }
        }
        for (; j < spline.bcoef.extent(1); ++j) {
            for (int i(0); i < spline.bcoef.extent(0); ++i) {
                spline.bcoef(i, j) = t_storage(j, i);
            }
        }
    }
    if (xmin_bc[1] == BoundCond::PERIODIC) {
        for (int i(0); i < spline.bcoef.extent(0); ++i) {
            for (int j(0); j < bspl[1]->degree; ++j) {
                spline.bcoef(i, j + bspl[1]->nbasis) = spline.bcoef(i, j);
            }
        }
    }
    if (xmin_bc[0] == BoundCond::PERIODIC) {
        for (int i(0); i < bspl[0]->degree; ++i) {
            for (int j(0); j < spline.bcoef.extent(1); ++j) {
                spline.bcoef(i + bspl[0]->nbasis, j) = spline.bcoef(i, j);
            }
        }
    }
}
