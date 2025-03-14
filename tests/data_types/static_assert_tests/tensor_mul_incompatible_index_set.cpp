// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include "../geometry_tensor.hpp"

#include "indexed_tensor.hpp"
#include "tensor.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

using Tensor2D_A = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;
using Tensor2D_B = Tensor<int, VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>>;
using Tensor2D_C = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;

Tensor2D_C compile_tensor_test_mul(Tensor2D_A const& A, Tensor2D_B const& B)
{
    return tensor_mul(index<'i', 'j'>(A), index<'j', 'k'>(B));
}
