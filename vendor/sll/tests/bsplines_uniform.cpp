#include <iosfwd>
#include <utility>
#include <vector>

#include <experimental/mdspan>

#include <ddc/MCoord>
#include <ddc/MDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <sll/bsplines_uniform.h>
#include <sll/deprecated/bsplines_uniform.h>
#include <sll/view.h>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using MeshX = UniformMesh<DimX>;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<MeshX>;

class BSplinesUniformTest : public ::testing::Test
{
protected:
    static constexpr std::size_t spline_degree = 2;
    static constexpr std::size_t npoints = 101;
    static constexpr RCoordX xmin = 0.;
    static constexpr RCoordX xmax = 2.;
    BSplines<MeshX, spline_degree> const bsplines {xmin, xmax, npoints};
    deprecated::UniformBSplines const
            old_bsplines {spline_degree, DimX::PERIODIC, xmin, xmax, npoints - 1};
};

TEST_F(BSplinesUniformTest, constructor)
{
    EXPECT_EQ(bsplines.degree(), spline_degree);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.rmin(), xmin);
    EXPECT_EQ(bsplines.rmax(), xmax);
    EXPECT_EQ(bsplines.npoints(), npoints);
    EXPECT_EQ(bsplines.ncells(), npoints - 1);
}

TEST_F(BSplinesUniformTest, comparison)
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

    bsplines.eval_basis(x, values, jmin);
    old_bsplines.eval_basis(x, old_values, old_jmin);
    EXPECT_EQ(jmin, old_jmin);
    for (std::size_t i = 0; i < values.extent(0); ++i) {
        EXPECT_EQ(values(i), old_values(i));
    }

    bsplines.eval_deriv(x, values, jmin);
    old_bsplines.eval_deriv(x, old_values, old_jmin);
    EXPECT_EQ(jmin, old_jmin);
    for (std::size_t i = 0; i < values.extent(0); ++i) {
        EXPECT_EQ(values(i), old_values(i));
    }
}
