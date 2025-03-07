// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "../geometry_tensor.hpp"

#include "indexed_tensor.hpp"
#include "tensor.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

TEST(TensorTest, Mul)
{
    using Tensor2D_A = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;
    using Tensor2D_B
            = Tensor<int, VectorIndexSet<R_cov, Theta_cov>, VectorIndexSet<R_cov, Theta_cov>>;
    using Tensor2D_C = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;
    Tensor2D_A A;
    Tensor2D_B B;
    ddcHelper::get<R, R>(A) = 1;
    ddcHelper::get<R, Theta>(A) = 0;
    ddcHelper::get<Theta, R>(A) = 2;
    ddcHelper::get<Theta, Theta>(A) = 4;
    ddcHelper::get<R, R_cov>(B) = 6;
    ddcHelper::get<R, Theta_cov>(B) = 8;
    ddcHelper::get<Theta, R_cov>(B) = 4;
    ddcHelper::get<Theta, Theta_cov>(B) = 5;
    Tensor2D_C C = tensor_mul(index<'i', 'j'>(A), index<'j', 'k'>(B));
}
