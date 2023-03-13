#include <iosfwd>
#include <utility>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>

#include <gtest/gtest.h>

namespace {

struct DimX
{
    [[maybe_unused]] static constexpr bool PERIODIC = true;
};
struct DimY
{
    [[maybe_unused]] static constexpr bool PERIODIC = false;
};

using CoordX = ddc::Coordinate<DimX>;
using IDimX = ddc::UniformPointSampling<DimX>;
using IndexX = ddc::DiscreteElement<IDimX>;
using BSplinesX = UniformBSplines<DimX, 3>;

using RCoordY = ddc::Coordinate<DimY>;
using MeshY = ddc::NonUniformPointSampling<DimY>;
using MCoordY = ddc::DiscreteElement<MeshY>;
using BSplinesY = NonUniformBSplines<DimY, 4>;

constexpr std::size_t ncells = 100;
constexpr CoordX xmin(0.);
constexpr CoordX xmax(2.);

} // namespace

TEST(ChunkBSplinesTest, Constructor)
{
    ddc::init_discrete_space<BSplinesX>(xmin, xmax, ncells);

    ddc::init_discrete_space<BSplinesY>(
            std::initializer_list<RCoordY> {RCoordY(0.1), RCoordY(0.4), RCoordY(1.0)});

    ddc::DiscreteElement<BSplinesX, BSplinesY> start(0, 0);
    ddc::DiscreteVector<BSplinesX, BSplinesY> size(ncells, ncells);
    ddc::DiscreteDomain<BSplinesX, BSplinesY> dom(start, size);

    ddc::Chunk<double, ddc::DiscreteDomain<BSplinesX, BSplinesY>> chunk(dom);
    auto view = chunk.span_view();

    for (ddc::DiscreteElement<BSplinesX> ibsx : ddc::get_domain<BSplinesX>(chunk)) {
        for (ddc::DiscreteElement<BSplinesY> ibsy : ddc::get_domain<BSplinesY>(chunk)) {
            view(ibsx, ibsy) = 1.0;
        }
    }
}
