#include <cstdint>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "deprecated/bsplines_non_uniform.h"

#include "bsplines_non_uniform.h"
#include "mdomain.h"
#include "non_uniform_mesh.h"

struct DimX
{
    static constexpr bool PERIODIC = true;
};
using NonUniformMeshX = NonUniformMesh<DimX>;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<DimX>;

TEST(BSplinesNonUniform, Constructor)
{
    std::vector<double> const breaks({0.0, 0.5, 1.0, 1.5, 2.0});
    NonUniformMeshX const mesh(breaks, MCoordX(0));
    MDomainImpl<NonUniformMeshX> const dom(mesh, MCoordX(5));

    std::integral_constant<std::size_t, 2> constexpr spline_degree;
    auto&& bsplines = make_bsplines(dom, spline_degree);

    EXPECT_EQ(bsplines.degree(), spline_degree.value);
    EXPECT_EQ(bsplines.is_periodic(), DimX::PERIODIC);
    EXPECT_EQ(bsplines.rmin(), 0.);
    EXPECT_EQ(bsplines.rmax(), 2.);
    EXPECT_EQ(bsplines.npoints(), 5);
    EXPECT_EQ(bsplines.ncells(), 4);
}

TEST(BSplinesNonUniform, Comparison)
{
    std::vector<double> const breaks({0.0, 0.5, 1.0, 1.5, 2.0});
    NonUniformMeshX const mesh(breaks, MCoordX(0));
    MDomainImpl<NonUniformMeshX> const dom(mesh, MCoordX(5));

    std::integral_constant<std::size_t, 2> constexpr spline_degree;
    auto&& bsplines = make_bsplines(dom, spline_degree);

    deprecated::NonUniformBSplines old_bsplines(spline_degree.value, dom);

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
