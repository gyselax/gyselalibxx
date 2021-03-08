#ifndef SPLINE_INTERPOLATORS_1D_H
#define SPLINE_INTERPOLATORS_1D_H
#include "bsplines.h"
#include "bsplines_non_uniform.h"
#include "bsplines_uniform.h"
#include "spline_1d.h"
#include <boundary_conditions.h>
#include <matrix.h>

class Spline_interpolator_1D {
    public:
        Spline_interpolator_1D(BSplines* bspl, BoundaryCondition xmin_bc, BoundaryCondition xmax_bc);
        ~Spline_interpolator_1D();
        const mdspan_1d& get_interp_points() const;
        void compute_interpolant(Spline_1D& spline,
                const mdspan_1d& vals,
                const mdspan_1d* derivs_xmin = nullptr,
                const mdspan_1d* derivs_xmax = nullptr) const;
        static int compute_num_cells(int degree, BoundaryCondition xmin, BoundaryCondition xmax, int nipts);
        const BoundaryCondition xmin_bc;
        const BoundaryCondition xmax_bc;
        const int nbc_xmin;
        const int nbc_xmax;
    private:
        void compute_interpolation_points_uniform( );
        void compute_interpolation_points_non_uniform( );
        void compute_block_sizes_uniform(int& kl, int& ku) const;
        void compute_block_sizes_non_uniform(int& kl, int& ku) const;

        void constructor_sanity_checks();
        void allocate_matrix(int kl, int ku);
        void compute_interpolant_degree1(Spline_1D& spline, const mdspan_1d& vals) const;
        void build_matrix_system();

        const BSplines* const bspl;
        const bool odd;
        const int offset;
        const double dx; // average cell size for normalization of derivatives
        double* interp_pts_ptr;
        mdspan_1d interp_pts;
        Matrix* matrix;
        static std::array<BoundaryCondition, 3> allowed_bcs;
};

extern "C"
{
Spline_interpolator_1D* new_spline_interpolator_1d(BSplines* bspl, BoundaryCondition xmin_bc, BoundaryCondition xmax_bc);
void free_spline_interpolator_1d(Spline_interpolator_1D* spl_interp);
int compute_num_cells(int degree, BoundaryCondition xmin, BoundaryCondition xmax, int nipts);
void compute_interpolant(const Spline_interpolator_1D* spl_interp, Spline_1D* spline,
        double* vals_ptr, int nvals, double* derivs_xmin_ptr = nullptr,
        int nderivs_xmin = 0, double* derivs_xmax_ptr = nullptr, int nderivs_xmax = 0);
void get_interp_points(Spline_interpolator_1D* spl_interp, double* interp_points, int npts);
int get_n_interp_points(Spline_interpolator_1D* spl_interp);
}
#endif // SPLINE_INTERPOLATORS_1D_H
