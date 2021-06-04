#pragma once

#include <array>
#include <memory>

#include "spline_2d.h"
#include "spline_builder_1d.h"

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

class SplineBuilder2D
{
private:
    const std::array<std::unique_ptr<const BSplines>, 2> bspl;
    std::array<SplineBuilder1D, 2> interp_1d;
    // TODO: Improve
    std::array<Spline1D, 2> spline_1d;

public:
    const std::array<BoundCond, 2> m_xmin_bc;
    const std::array<BoundCond, 2> m_xmax_bc;
    const std::array<int, 2> m_nbc_xmin;
    const std::array<int, 2> m_nbc_xmax;

public:
    SplineBuilder2D() = delete;
    SplineBuilder2D(
            std::array<std::unique_ptr<const BSplines>, 2> bspl,
            std::array<BoundCond, 2> xmin_bc,
            std::array<BoundCond, 2> xmax_bc);
    SplineBuilder2D(const SplineBuilder2D& x) = delete;
    SplineBuilder2D(SplineBuilder2D&& x) = delete;
    ~SplineBuilder2D() = default;
    SplineBuilder2D& operator=(const SplineBuilder2D& x) = delete;
    SplineBuilder2D& operator=(SplineBuilder2D&& x) = delete;

    std::array<const DSpan1D, 2> get_interp_points() const;

    void compute_interpolant(
            Spline2D const& spline,
            DSpan2D const& vals,
            Boundary_data_2d boundary_data) const;

    void compute_interpolant(Spline2D const& spline, DSpan2D const& vals) const;

    static std::array<int, 2> compute_num_cells(
            std::array<int, 2> degree,
            std::array<BoundCond, 2> xmin,
            std::array<BoundCond, 2> xmax,
            std::array<int, 2> nipts);

    std::array<BoundCond, 2> const& xmin_bc() const noexcept
    {
        return m_xmin_bc;
    }

    std::array<BoundCond, 2> const& xmax_bc() const noexcept
    {
        return m_xmax_bc;
    }

    std::array<int, 2> const& nbc_xmin() const noexcept
    {
        return m_nbc_xmin;
    }

    std::array<int, 2> const& nbc_xmax() const noexcept
    {
        return m_nbc_xmax;
    }

private:
    void compute_interpolant_boundary_done(Spline2D const& spline, DSpan2D const& vals) const;
};