#ifndef SPLINE_INTERPOLATORS_1D_H
#define SPLINE_INTERPOLATORS_1D_H
#include <boundary_conditions.h>
#include <matrix.h>

#include <memory>

#include "bsplines.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"
#include "spline_1d.h"

class Spline_interpolator_1D {
public:
    Spline_interpolator_1D(const BSplines& bspl,
        BoundaryCondition xmin_bc,
        BoundaryCondition xmax_bc);
    const mdspan_1d& get_interp_points() const;
    void compute_interpolant(Spline_1D& spline,
        const mdspan_1d& vals,
        const mdspan_1d* derivs_xmin = nullptr,
        const mdspan_1d* derivs_xmax = nullptr) const;
    static int compute_num_cells(int degree,
        BoundaryCondition xmin,
        BoundaryCondition xmax,
        int nipts);
    const BoundaryCondition xmin_bc;
    const BoundaryCondition xmax_bc;
    const int nbc_xmin;
    const int nbc_xmax;

private:
    void compute_interpolation_points_uniform();
    void compute_interpolation_points_non_uniform();
    void compute_block_sizes_uniform(int& kl, int& ku) const;
    void compute_block_sizes_non_uniform(int& kl, int& ku) const;

    void constructor_sanity_checks();
    void allocate_matrix(int kl, int ku);
    void compute_interpolant_degree1(Spline_1D& spline,
        const mdspan_1d& vals) const;
    void build_matrix_system();

    const BSplines& bspl;
    const bool odd;
    const int offset;
    const double dx; // average cell size for normalization of derivatives
    std::unique_ptr<double[]> interp_pts_ptr;
    mdspan_1d interp_pts;
    std::unique_ptr<Matrix> matrix;
    static std::array<BoundaryCondition, 3> allowed_bcs;
};

#endif // SPLINE_INTERPOLATORS_1D_H
