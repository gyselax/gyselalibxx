// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "../geometry_tensor.hpp"

#include "indexed_tensor.hpp"
#include "tensor.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

TEST(TensorTest, Index)
{
    using Tensor2D_A = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;
    Tensor2D_A A;
    ddcHelper::get<R, X>(A) = 1;
    ddcHelper::get<R, Y>(A) = 0;
    ddcHelper::get<Theta, X>(A) = 2;
    ddcHelper::get<Theta, Y>(A) = 4;
}
