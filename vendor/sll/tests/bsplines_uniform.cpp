#include <iosfwd>
#include <utility>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_uniform.hpp>
#include <sll/deprecated/bsplines_uniform.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using IDimX = UniformPointSampling<DimX>;
using CoordX = Coordinate<DimX>;
using IndexX = DiscreteElement<IDimX>;

class BSplinesUniformTest : public ::testing::Test
{
protected:
    static constexpr std::size_t spline_degree = 2;
    static constexpr std::size_t ncells = 100;
    static constexpr CoordX xmin = CoordX(0.);
    static constexpr CoordX xmax = CoordX(2.);
    UniformBSplines<DimX, spline_degree> const bsplines {xmin, xmax, ncells};
    deprecated::UniformBSplines const
            old_bsplines {spline_degree, DimX::PERIODIC, xmin, xmax, ncells};
};

TEST_F(BSplinesUniformTest, Constructor)
{
    EXPECT_EQ(bsplines.degree(), spline_degree);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.rmin(), xmin);
    EXPECT_EQ(bsplines.rmax(), xmax);
    EXPECT_EQ(bsplines.ncells(), ncells);
}

TEST_F(BSplinesUniformTest, Comparison)
{
    EXPECT_EQ(bsplines.degree(), old_bsplines.degree());
    EXPECT_EQ(bsplines.is_radial(), old_bsplines.radial());
    EXPECT_EQ(bsplines.is_periodic(), old_bsplines.is_periodic());
    EXPECT_EQ(bsplines.is_uniform(), old_bsplines.is_uniform());
    EXPECT_EQ(bsplines.nbasis(), old_bsplines.nbasis());
    EXPECT_EQ(bsplines.ncells(), old_bsplines.ncells());
    EXPECT_EQ(bsplines.rmin(), old_bsplines.xmin());
    EXPECT_EQ(bsplines.rmax(), old_bsplines.xmax());
    EXPECT_EQ(bsplines.length(), old_bsplines.length());

    std::vector<double> values_data(bsplines.degree() + 1);
    DSpan1D values(values_data.data(), values_data.size());
    int jmin;
    std::vector<double> old_values_data(old_bsplines.degree() + 1);
    DSpan1D old_values(old_values_data.data(), old_values_data.size());
    int old_jmin;

    double const x = 1.07;

    bsplines.eval_basis(values, jmin, x);
    old_bsplines.eval_basis(x, old_values, old_jmin);
    EXPECT_EQ(jmin, old_jmin);
    for (std::size_t i = 0; i < values.extent(0); ++i) {
        EXPECT_EQ(values(i), old_values(i));
    }

    bsplines.eval_deriv(values, jmin, x);
    old_bsplines.eval_deriv(x, old_values, old_jmin);
    EXPECT_EQ(jmin, old_jmin);
    for (std::size_t i = 0; i < values.extent(0); ++i) {
        EXPECT_EQ(values(i), old_values(i));
    }
}
