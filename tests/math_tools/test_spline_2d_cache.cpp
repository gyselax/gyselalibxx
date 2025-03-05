// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "mesh_builder.hpp"
#include "spline_builder_2d_cache.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};

struct Y
{
    static bool constexpr PERIODIC = false;
};

using CoordX = Coord<X>;
using CoordY = Coord<Y>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
struct BSplinesY : ddc::UniformBSplines<Y, 3>
{
};

auto constexpr SplineXBoundary = ddc::BoundCond::GREVILLE;
auto constexpr SplineYBoundary = ddc::BoundCond::GREVILLE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;

struct GridX : NonUniformGridBase<X>
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};

using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

using IdxRangeXY = IdxRange<GridX, GridY>;

using SplineXYBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        SplineXBoundary,
        SplineXBoundary,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;

using IdxRangeBSXY = SplineXYBuilder::batched_spline_domain_type;

// Dummy spline builder class for keeping track of how many times builder is called by spline 2d cache
class DummySplineBuilder2D
{
public:
    using continuous_dimension_type1 = X;
    using continuous_dimension_type2 = Y;
    using batched_interpolation_domain_type = ddc::DiscreteDomain<GridX, GridY>;
    using batched_spline_domain_type = ddc::DiscreteDomain<BSplinesX, BSplinesY>;
    using DFieldSplineCoeff = DField<batched_spline_domain_type>;
    using DConstFieldCoeff = DConstField<batched_interpolation_domain_type>;

    IdxRangeBSXY m_spline_domain;
    int mutable m_builder_call_counter;

public:
    /**
  * @brief Creates an instance of the class DummySplineBuilder2D.
  * This class contains all necessary elements to be passed as a 
  * template parameter type to the SplineBuilder2DCache class.
  *
  * @param spline_domain The domain on which splines are constructed.  
  */
    DummySplineBuilder2D(IdxRangeBSXY spline_domain)
        : m_spline_domain(spline_domain)
        , m_builder_call_counter(0)
    {
    }

    /**
     * @brief Get the whole domain on which spline coefficients are defined.
     *
     * @return The domain for the spline coefficients.
     */
    batched_spline_domain_type batched_spline_domain() const noexcept
    {
        return m_spline_domain;
    }

    /**
     * @brief increments the counter that tracks the number of calls to the builder. 
     */
    void operator()(DFieldSplineCoeff, DConstFieldCoeff) const
    {
        m_builder_call_counter++;
    }
};

TEST(SplineBuilder2DCache, CountBuilderCalls)
{
    int n_elems_x(10);
    int n_elems_y(20);

    Coord<X> const x_min(0.0);
    Coord<X> const x_max(1.0);
    IdxStepX x_ncells(n_elems_x);

    Coord<Y> const y_min(0.0);
    Coord<Y> const y_max(2.0);
    IdxStepY y_ncells(n_elems_y);

    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX idxrange_x(SplineInterpPointsX::get_domain<GridX>());

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_ncells);
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    IdxRangeY idxrange_y(SplineInterpPointsY::get_domain<GridY>());

    IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);

    IdxRangeBSXY idxrange_bs_xy(
            ddc::discrete_space<BSplinesX>().full_domain(),
            ddc::discrete_space<BSplinesY>().full_domain());

    DummySplineBuilder2D dummy_builder(idxrange_bs_xy);
    SplineBuilder2DCache<DummySplineBuilder2D> builder_cache(dummy_builder);

    DFieldMem<IdxRangeXY> field_xy_alloc(idxrange_xy);
    DField<IdxRangeXY> field_xy = get_field(field_xy_alloc);
    ddc::parallel_fill(field_xy, 1.);

    builder_cache.compute_coeffs<X>(get_const_field(field_xy));
    builder_cache.compute_coeffs<Y>(get_const_field(field_xy));
    EXPECT_EQ(dummy_builder.m_builder_call_counter, 1);

    builder_cache.compute_coeffs<Y>(get_const_field(field_xy));
    builder_cache.compute_coeffs<X>(get_const_field(field_xy));
    EXPECT_EQ(dummy_builder.m_builder_call_counter, 2);
}

} // namespace
