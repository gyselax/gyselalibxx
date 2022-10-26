#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/null_boundary_value.hpp>
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

using IDimX = SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC>::
        interpolation_mesh_type;

using IndexX = DiscreteElement<IDimX>;
using DVectX = DiscreteVector<IDimX>;
using BsplIndexX = DiscreteElement<BSplinesX>;
using SplineX = Chunk<double, DiscreteDomain<BSplinesX>>;
using FieldX = Chunk<double, DiscreteDomain<IDimX>>;
using CoordX = Coordinate<DimX>;

TEST(PeriodicSplineBuilderOrderTest, OrderedPoints)
{
    CoordX constexpr x0(0.);
    CoordX constexpr xN(1.);
    std::size_t constexpr ncells = 10;

    // 1. Create BSplines
    int constexpr npoints(ncells + 1);
    std::vector<double> d_breaks({0, 0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
    std::vector<CoordX> breaks(npoints);
    for (std::size_t i(0); i < npoints; ++i) {
        breaks[i] = CoordX(d_breaks[i]);
    }
    init_discrete_space<BSplinesX>(breaks);

    // 2. Create a SplineBuilder over BSplines using some boundary conditions
    SplineBuilder<BSplinesX, BoundCond::PERIODIC, BoundCond::PERIODIC> spline_builder;
    auto interpolation_domain = spline_builder.interpolation_domain();

    double last(coordinate(interpolation_domain.front()));
    double current;
    for (IndexX const ix : interpolation_domain) {
        current = coordinate(ix);
        ASSERT_LE(current, discrete_space<BSplinesX>().rmax());
        ASSERT_GE(current, discrete_space<BSplinesX>().rmin());
        ASSERT_LE(last, current);
        last = current;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::ScopeGuard scope(argc, argv);
    return RUN_ALL_TESTS();
}
