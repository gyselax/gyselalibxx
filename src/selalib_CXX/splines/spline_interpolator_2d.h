#ifndef SPLINE_INTERPOLATORS_2D_H
#define SPLINE_INTERPOLATORS_2D_H
#include "spline_2d.h"
#include "spline_interpolator_1d.h"
#include <array>

struct Boundary_data_2d
{
    mdspan_2d* derivs_x1_min = nullptr;
    mdspan_2d* derivs_x1_max = nullptr;
    mdspan_2d* derivs_x2_min = nullptr;
    mdspan_2d* derivs_x2_max = nullptr;
    mdspan_2d* mixed_derivs_a = nullptr;
    mdspan_2d* mixed_derivs_b = nullptr;
    mdspan_2d* mixed_derivs_c = nullptr;
    mdspan_2d* mixed_derivs_d = nullptr;
};

class Spline_interpolator_2D {
    public:
        Spline_interpolator_2D(std::array<BSplines*,2> bspl,
                std::array<BoundaryCondition,2> xmin_bc,
                std::array<BoundaryCondition,2> xmax_bc);

        ~Spline_interpolator_2D() {}

        std::array<const mdspan_1d,2> get_interp_points() const;

        void compute_interpolant(Spline_2D const& spline,
                mdspan_2d const& vals,
                Boundary_data_2d boundary_data) const;

        void compute_interpolant(Spline_2D const& spline,
                mdspan_2d const& vals) const;

        static std::array<int,2> compute_num_cells(std::array<int,2> degree,
                std::array<BoundaryCondition,2> xmin,
                std::array<BoundaryCondition,2> xmax,
                std::array<int,2> nipts);
    private:
        void compute_interpolant_boundary_done(Spline_2D const& spline,
                mdspan_2d const& vals) const;

        const std::array<BSplines*, 2> bspl;
        std::array<Spline_interpolator_1D,2> interp_1d;
        // TODO: Improve
        std::array<Spline_1D,2> spline_1d;
    public:
        const std::array<BoundaryCondition, 2> xmin_bc;
        const std::array<BoundaryCondition, 2> xmax_bc;
        const std::array<int, 2> nbc_xmin;
        const std::array<int, 2> nbc_xmax;
};

extern "C"
{
Spline_interpolator_2D* new_spline_interpolator_2d(BSplines* bspl_1,
        BoundaryCondition xmin_bc_1, BoundaryCondition xmax_bc_1,
        BSplines* bspl_2, BoundaryCondition xmin_bc_2, BoundaryCondition xmax_bc_2);
void free_spline_interpolator_2d(Spline_interpolator_2D* spl_interp);
void compute_num_cells_2d(int degree_1, BoundaryCondition xmin_1, BoundaryCondition xmax_1, int nipts_1,
        int degree_2, BoundaryCondition xmin_2, BoundaryCondition xmax_2, int nipts_2,
        int* ncell_1, int* ncell_2);
void compute_interpolant_2d(const Spline_interpolator_2D* spl_interp, Spline_2D* spline,
        double* vals_ptr, int nvals_1, int nvals_2,
        double* derivs_x1_min,  int n_derivs_x1_min,  int n_derivs_x1_min_2,
        double* derivs_x1_max,  int n_derivs_x1_max,  int n_derivs_x1_max_2,
        double* derivs_x2_min,  int n_derivs_x2_min,  int n_derivs_x2_min_2,
        double* derivs_x2_max,  int n_derivs_x2_max,  int n_derivs_x2_max_2,
        double* mixed_derivs_a, int n_mixed_derivs_a, int n_mixed_derivs_a_2,
        double* mixed_derivs_b, int n_mixed_derivs_b, int n_mixed_derivs_b_2,
        double* mixed_derivs_c, int n_mixed_derivs_c, int n_mixed_derivs_c_2,
        double* mixed_derivs_d, int n_mixed_derivs_d, int n_mixed_derivs_d_2);
void get_interp_points_2d(Spline_interpolator_2D* spl_interp, double* interp_points_1, int npts_1,
        double* interp_points_2, int npts_2);
void get_n_interp_points_2d(Spline_interpolator_2D* spl_interp, int* npts_1, int* npts_2);
}
#endif // SPLINE_INTERPOLATORS_2D_H
