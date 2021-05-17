#ifndef SPLINE_INTERPOLATORS_1D_H
#define SPLINE_INTERPOLATORS_1D_H
#include <memory>

#include <boundary_conditions.h>
#include <matrix.h>

#include "bsplines.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"
#include "spline_1d.h"

class Spline_interpolator_1D
{
public:
    Spline_interpolator_1D(const BSplines& bspl, BoundCond xmin_bc, BoundCond xmax_bc);

    const DSpan1D& get_interp_points() const;

    void compute_interpolant(
            Spline_1D& spline,
            const DSpan1D& vals,
            const DSpan1D* derivs_xmin = nullptr,
            const DSpan1D* derivs_xmax = nullptr) const;

    static int compute_num_cells(int degree, BoundCond xmin, BoundCond xmax, int nipts);

    const BoundCond xmin_bc;

    const BoundCond xmax_bc;

    const int nbc_xmin;

    const int nbc_xmax;

private:
    void compute_interpolation_points_uniform();
    void compute_interpolation_points_non_uniform();
    void compute_block_sizes_uniform(int& kl, int& ku) const;
    void compute_block_sizes_non_uniform(int& kl, int& ku) const;

    void constructor_sanity_checks();
    void allocate_matrix(int kl, int ku);
    void compute_interpolant_degree1(Spline_1D& spline, const DSpan1D& vals) const;
    void build_matrix_system();

    const BSplines& bspl;

    // bspline stuff: TODO move
    const bool odd; // bspl.degree % 2 == 1
    const int offset; // bspl.periodic ? bspl.degree / 2 : 0
    const double dx; // average cell size for normalization of derivatives

    // mesh info: TODO use Mesh
    std::unique_ptr<double[]> interp_pts_ptr;
    DSpan1D interp_pts;

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

    // hand-made inheritance
    static std::array<BoundCond, 3> allowed_bcs;
};

#endif // SPLINE_INTERPOLATORS_1D_H
