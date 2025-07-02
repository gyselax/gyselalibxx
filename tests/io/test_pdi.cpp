#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "pdi_helper.hpp"

namespace {

struct GridX
{
};

struct Density;
struct Velocity;

using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  moment_density_extents : { type: array, subtype: int64, size: 1 }
  moment_velocity_extents : { type: array, subtype: int64, size: 1 }
  moment_density:
    type: array
    subtype: double
    size: [ '$moment_density_extents[0]' ]
  moment_velocity:
    type: array
    subtype: double
    size: [ '$moment_velocity_extents[0]' ]

data:

plugins:
  decl_hdf5:
    - file: 'test_pdi.h5'
      on_event: [write_moments]
      collision_policy: replace_and_warn
      write: [moment_density, moment_velocity, moment_density_extents, moment_velocity_extents]
    - file: 'test_pdi.h5'
      on_event: [read_moments_extents]
      read: [moment_density_extents, moment_velocity_extents]
    - file: 'test_pdi.h5'
      on_event: [read_moments]
      read: [moment_density, moment_velocity]
  #trace: ~
)PDI_CFG";

TEST(PDI, VectorFieldReadWrite)
{
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);
    IdxRangeX idxrange(IdxX(0), IdxStepX(10));
    host_t<DVectorFieldMem<IdxRangeX, VectorIndexSet<Density, Velocity>>> my_moments(idxrange);

    ddc::parallel_fill(ddcHelper::get<Density>(my_moments), 1.0);
    ddc::parallel_fill(ddcHelper::get<Velocity>(my_moments), -9.5);

    PDI_expose_vector_field("moment_", get_const_field(my_moments), "density", "velocity");
    ddc::PdiEvent("write_moments");

    std::vector<double> density;
    std::vector<double> velocity;

    PDI_get_arrays("read_moments", "moment_density", density, "moment_velocity", velocity);

    EXPECT_EQ(density.size(), idxrange.size());
    EXPECT_EQ(velocity.size(), idxrange.size());

    for (int i(0); i < idxrange.size(); ++i) {
        EXPECT_DOUBLE_EQ(density[i], ddcHelper::get<Density>(my_moments)(IdxStepX(i)));
        EXPECT_DOUBLE_EQ(velocity[i], ddcHelper::get<Velocity>(my_moments)(IdxStepX(i)));
    }

    PDI_finalize();
}

} // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    return RUN_ALL_TESTS();
}
