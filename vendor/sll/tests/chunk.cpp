#include <iosfwd>
#include <utility>
#include <vector>

#include <experimental/mdspan>

#include <ddc/Chunk>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/UniformDiscretization>

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
using IDimX = UniformDiscretization<DimX>;
using IndexX = DiscreteCoordinate<IDimX>;
using BSplinesX = UniformBSplines<DimX, 3>;

using RCoordY = Coordinate<DimY>;
using MeshY = NonUniformDiscretization<DimY>;
using MCoordY = DiscreteCoordinate<MeshY>;
using BSplinesY = NonUniformBSplines<DimY, 4>;

constexpr std::size_t ncells = 100;
constexpr CoordX xmin(0.);
constexpr CoordX xmax(2.);
BSplinesX const bsplinesx {xmin, xmax, ncells};

BSplinesY const bsplinesy {RCoordY(0.1), RCoordY(0.4), RCoordY(1.0)};

} // namespace

TEST(ChunkBSplinesTest, Constructor)
{
    DiscreteVector<BSplinesX, BSplinesY> size(ncells, ncells);
    DiscreteDomain<BSplinesX, BSplinesY> dom(size);

    Chunk<double, DiscreteDomain<BSplinesX, BSplinesY>> chunk(dom);
    auto view = chunk.span_view();

    for (DiscreteCoordinate<BSplinesX> ibsx : get_domain<BSplinesX>(chunk)) {
        for (DiscreteCoordinate<BSplinesY> ibsy : get_domain<BSplinesY>(chunk)) {
            view(ibsx, ibsy) = 1.0;
        }
    }
}
