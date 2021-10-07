#include <iosfwd>
#include <utility>
#include <vector>

#include <experimental/mdspan>

#include <ddc/MCoord>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>
#include <ddc/block.hpp>

#include <sll/bsplines_non_uniform.h>
#include <sll/bsplines_uniform.h>

#include <gtest/gtest.h>

struct DimX
{
    static constexpr bool PERIODIC = true;
};
struct DimY
{
    static constexpr bool PERIODIC = false;
};

using MeshX = UniformMesh<DimX>;
using MeshY = NonUniformMesh<DimY>;
using RCoordX = RCoord<DimX>;
using MCoordX = MCoord<MeshX>;
using RCoordY = RCoord<DimY>;
using MCoordY = MCoord<MeshY>;

class BlockBSplinesUniformTest : public ::testing::Test
{
protected:
    static constexpr std::size_t ncells = 100;
    static constexpr RCoordX xmin = RCoordX(0.);
    static constexpr RCoordX xmax = RCoordX(2.);
    using BSplinesX = UniformBSplines<DimX, 2>;
    BSplinesX const bsplinesx {xmin, xmax, ncells};

    using BSplinesY = NonUniformBSplines<DimY, 4>;
    BSplinesY const bsplinesy {RCoordY(0.1), RCoordY(0.4), RCoordY(1.0)};
};

TEST_F(BlockBSplinesUniformTest, constructor)
{
    MLength<BSplinesX, BSplinesY> size(ncells, ncells);
    ProductMDomain<BSplinesX, BSplinesY> dom(bsplinesx, bsplinesy, size);

    Block<double, ProductMDomain<BSplinesX, BSplinesY>> block(dom);
    auto view = block.view();

    for (MCoord<BSplinesX> ibsx : get_domain<BSplinesX>(block)) {
        for (MCoord<BSplinesY> ibsy : get_domain<BSplinesY>(block)) {
            view(ibsx, ibsy) = 1.0;
        }
    }
}
