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

using CoordX = Coordinate<DimX>;
using IDimX = UniformPointSampling<DimX>;
using IndexX = DiscreteElement<IDimX>;
using BSplinesX = UniformBSplines<DimX, 3>;

using RCoordY = Coordinate<DimY>;
using MeshY = NonUniformPointSampling<DimY>;
using MCoordY = DiscreteElement<MeshY>;
using BSplinesY = NonUniformBSplines<DimY, 4>;

constexpr std::size_t ncells = 100;
constexpr CoordX xmin(0.);
constexpr CoordX xmax(2.);
BSplinesX::Impl<Kokkos::HostSpace> const bsplinesx {xmin, xmax, ncells};

BSplinesY::Impl<Kokkos::HostSpace> const bsplinesy {RCoordY(0.1), RCoordY(0.4), RCoordY(1.0)};

} // namespace

TEST(ChunkBSplinesTest, Constructor)
{
    DiscreteElement<BSplinesX, BSplinesY> start(0, 0);
    DiscreteVector<BSplinesX, BSplinesY> size(ncells, ncells);
    DiscreteDomain<BSplinesX, BSplinesY> dom(start, size);

    Chunk<double, DiscreteDomain<BSplinesX, BSplinesY>> chunk(dom);
    auto view = chunk.span_view();

    for (DiscreteElement<BSplinesX> ibsx : get_domain<BSplinesX>(chunk)) {
        for (DiscreteElement<BSplinesY> ibsy : get_domain<BSplinesY>(chunk)) {
            view(ibsx, ibsy) = 1.0;
        }
    }
}
