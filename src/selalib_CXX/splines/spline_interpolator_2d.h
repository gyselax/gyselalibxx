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
        Spline_interpolator_2D(std::array<std::unique_ptr<const BSplines>,2> bspl,
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

        const std::array<std::unique_ptr<const BSplines>, 2> bspl;
        std::array<Spline_interpolator_1D,2> interp_1d;
        // TODO: Improve
        std::array<Spline_1D,2> spline_1d;
    public:
        const std::array<BoundaryCondition, 2> xmin_bc;
        const std::array<BoundaryCondition, 2> xmax_bc;
        const std::array<int, 2> nbc_xmin;
        const std::array<int, 2> nbc_xmax;
};
#endif // SPLINE_INTERPOLATORS_2D_H
