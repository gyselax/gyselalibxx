#pragma once

#include <array>
#include <memory>

#include "sll/deprecated/boundary_conditions.h"
#include "sll/matrix.h"
#include "sll/view.h"

namespace deprecated {
class BSplines;
class Spline1D;

class SplineBuilder1D
{
private:
    // hand-made inheritance
    static std::array<BoundCond, 3> allowed_bcs;

private:
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

    const BoundCond m_xmin_bc;

    const BoundCond m_xmax_bc;

    const int m_nbc_xmin;

    const int m_nbc_xmax;

public:
    SplineBuilder1D() = delete;
    SplineBuilder1D(const BSplines& bspl, BoundCond xmin_bc, BoundCond xmax_bc);
    SplineBuilder1D(const SplineBuilder1D& x) = delete;
    SplineBuilder1D(SplineBuilder1D&& x) = delete;
    ~SplineBuilder1D() = default;
    SplineBuilder1D& operator=(const SplineBuilder1D& x) = delete;
    SplineBuilder1D& operator=(SplineBuilder1D&& x) = delete;

    const DSpan1D& get_interp_points() const;

    void compute_interpolant(
            Spline1D& spline,
            const DSpan1D& vals,
            const DSpan1D* derivs_xmin = nullptr,
            const DSpan1D* derivs_xmax = nullptr) const;

    static int compute_num_cells(int degree, BoundCond xmin, BoundCond xmax, int nipts);

    BoundCond xmin_bc() const noexcept
    {
        return m_xmin_bc;
    }

    BoundCond xmax_bc() const noexcept
    {
        return m_xmax_bc;
    }

    int nbc_xmin() const noexcept
    {
        return m_nbc_xmin;
    }

    int nbc_xmax() const noexcept
    {
        return m_nbc_xmax;
    }

private:
    void compute_interpolation_points_uniform();
    void compute_interpolation_points_non_uniform();
    void compute_block_sizes_uniform(int& kl, int& ku) const;
    void compute_block_sizes_non_uniform(int& kl, int& ku) const;

    void constructor_sanity_checks() const;
    void allocate_matrix(int kl, int ku);
    void compute_interpolant_degree1(Spline1D& spline, const DSpan1D& vals) const;
    void build_matrix_system();
};

} // namespace deprecated
