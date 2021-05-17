#ifndef SPLINE_INTERPOLATORS_2D_H
#define SPLINE_INTERPOLATORS_2D_H
#include <array>

#include "spline_2d.h"
#include "spline_interpolator_1d.h"

struct Boundary_data_2d
{
    DSpan2D* derivs_x1_min = nullptr;
    DSpan2D* derivs_x1_max = nullptr;
    DSpan2D* derivs_x2_min = nullptr;
    DSpan2D* derivs_x2_max = nullptr;
    DSpan2D* mixed_derivs_a = nullptr;
    DSpan2D* mixed_derivs_b = nullptr;
    DSpan2D* mixed_derivs_c = nullptr;
    DSpan2D* mixed_derivs_d = nullptr;
};

class Spline_interpolator_2D
{
public:
    Spline_interpolator_2D(
            std::array<std::unique_ptr<const BSplines>, 2> bspl,
            std::array<BoundCond, 2> xmin_bc,
            std::array<BoundCond, 2> xmax_bc);

    ~Spline_interpolator_2D() {}

    std::array<const DSpan1D, 2> get_interp_points() const;

    void compute_interpolant(
            Spline_2D const& spline,
            DSpan2D const& vals,
            Boundary_data_2d boundary_data) const;

    void compute_interpolant(Spline_2D const& spline, DSpan2D const& vals) const;

    static std::array<int, 2> compute_num_cells(
            std::array<int, 2> degree,
            std::array<BoundCond, 2> xmin,
            std::array<BoundCond, 2> xmax,
            std::array<int, 2> nipts);

private:
    void compute_interpolant_boundary_done(Spline_2D const& spline, DSpan2D const& vals) const;

    const std::array<std::unique_ptr<const BSplines>, 2> bspl;
    std::array<Spline_interpolator_1D, 2> interp_1d;
    // TODO: Improve
    std::array<Spline_1D, 2> spline_1d;

public:
    const std::array<BoundCond, 2> xmin_bc;
    const std::array<BoundCond, 2> xmax_bc;
    const std::array<int, 2> nbc_xmin;
    const std::array<int, 2> nbc_xmax;
};
#endif // SPLINE_INTERPOLATORS_2D_H
