#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/MCoord>
#include <ddc/NonUniformMesh>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>
#include <ddc/RCoord>
#include <ddc/TaggedVector>

#include <sll/bsplines_non_uniform.h>
#include <sll/deprecated/bsplines_non_uniform.h>
#include <sll/view.h>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using MeshX = NonUniformMesh<DimX>;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<MeshX>;

class BSplinesNonUniformTest : public ::testing::Test
{
protected:
    static constexpr std::size_t spline_degree = 2;
    std::vector<double> const breaks {0.0, 0.5, 1.0, 1.5, 2.0};
    MeshX const mesh_x = MeshX(breaks);
    ProductMesh<MeshX> mesh = ProductMesh<MeshX>(mesh_x);
    ProductMDomain<MeshX> const dom = ProductMDomain(mesh, MCoordX(mesh_x.size() - 1));
    BSplines<MeshX, spline_degree> const bsplines = BSplines<MeshX, 2>(dom);
};

TEST_F(BSplinesNonUniformTest, constructor)
{
    EXPECT_EQ(bsplines.degree(), spline_degree);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.rmin(), 0.);
    EXPECT_EQ(bsplines.rmax(), 2.);
    EXPECT_EQ(bsplines.npoints(), 5);
    EXPECT_EQ(bsplines.ncells(), 4);
}

TEST_F(BSplinesNonUniformTest, comparison)
{
    deprecated::NonUniformBSplines old_bsplines(spline_degree, DimX::PERIODIC, breaks);

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
