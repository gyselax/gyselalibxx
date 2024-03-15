#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_boundary_conditions.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

using BSplinesX = NonUniformBSplines<DimX, s_degree_x>;

using GrevillePoints
        = GrevilleInterpolationPoints<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using IDimX = GrevillePoints::interpolation_mesh_type;

using IndexX = ddc::DiscreteElement<IDimX>;
using DVectX = ddc::DiscreteVector<IDimX>;
using BsplIndexX = ddc::DiscreteElement<BSplinesX>;
using SplineX = ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX>>;
using FieldX = ddc::Chunk<double, ddc::DiscreteDomain<IDimX>>;
using CoordX = ddc::Coordinate<DimX>;

TEST(PeriodicSplineBuilderOrderTest, OrderedPoints)
{
    std::size_t constexpr ncells = 10;

    // 1. Create BSplines
    int constexpr npoints(ncells + 1);
    std::vector<double> d_breaks({0, 0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
    std::vector<CoordX> breaks(npoints);
    for (std::size_t i(0); i < npoints; ++i) {
        breaks[i] = CoordX(d_breaks[i]);
    }
    ddc::init_discrete_space<BSplinesX>(breaks);

    // 2. Create the interpolation domain
    ddc::init_discrete_space<IDimX>(GrevillePoints::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain(GrevillePoints::get_domain());

    double last(ddc::coordinate(interpolation_domain.front()));
    double current;
    for (IndexX const ix : interpolation_domain) {
        current = ddc::coordinate(ix);
        ASSERT_LE(current, ddc::discrete_space<BSplinesX>().rmax());
        ASSERT_GE(current, ddc::discrete_space<BSplinesX>().rmin());
        ASSERT_LE(last, current);
        last = current;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
