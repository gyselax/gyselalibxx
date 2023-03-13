#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/deprecated/bsplines_non_uniform.hpp>
#include <sll/view.hpp>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using IDimX = ddc::NonUniformPointSampling<DimX>;
using CoordX = ddc::Coordinate<DimX>;
using IndexX = ddc::DiscreteElement<IDimX>;

class BSplinesNonUniformTest : public ::testing::Test
{
protected:
    static constexpr std::size_t spline_degree = 2;
    std::vector<double> const breaks {0.0, 0.5, 1.0, 1.5, 2.0};
    NonUniformBSplines<DimX, spline_degree> const bsplines {breaks};
    deprecated::NonUniformBSplines old_bsplines {spline_degree, DimX::PERIODIC, breaks};
};

TEST_F(BSplinesNonUniformTest, Constructor)
{
    EXPECT_EQ(bsplines.degree(), spline_degree);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.rmin(), 0.);
    EXPECT_EQ(bsplines.rmax(), 2.);
    EXPECT_EQ(bsplines.npoints(), 5);
    EXPECT_EQ(bsplines.ncells(), 4);
}

TEST_F(BSplinesNonUniformTest, Comparison)
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
