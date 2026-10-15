#include <gtest/gtest.h>

#include "matrix_utils.hpp"

namespace {

void fill_strided_values(
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> data,
        int const stride)
{
    Kokkos::parallel_for(
            "fill",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, data.extent(0)),
            KOKKOS_LAMBDA(int const i) {
                // Only every stride-th value is part of the vector, the others are larger
                // to check that they are ignored.
                data(i) = (i % stride == 0) ? -(i / stride) - 1.5 : 100.;
            });
}

} // namespace

TEST(MatrixUtils, InfNorm)
{
    int const n = 5;
    int const stride = 3;
    std::shared_ptr const gko_exec
            = gko::ext::kokkos::create_executor(Kokkos::DefaultExecutionSpace());
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            data("data", n * stride);
    fill_strided_values(data, stride);
    auto vec = gko::matrix::Dense<double>::
            create(gko_exec,
                   gko::dim<2>(n, 1),
                   gko::array<double>::view(gko_exec, data.span(), data.data()),
                   stride);
    EXPECT_DOUBLE_EQ(inf_norm(vec.get()), n + 0.5);
}
