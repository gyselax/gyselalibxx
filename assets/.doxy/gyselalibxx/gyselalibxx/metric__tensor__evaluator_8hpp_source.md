

# File metric\_tensor\_evaluator.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**metric\_tensor\_evaluator.hpp**](metric__tensor__evaluator_8hpp.md)

[Go to the documentation of this file](metric__tensor__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

template <class Mapping, class PositionCoordinate = typename Mapping::CoordArg>
class MetricTensorEvaluator
{
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_jacobian_v<Mapping, PositionCoordinate>);
    static_assert(
            (is_covariant_vector_index_set_v<ddc::to_type_seq_t<typename Mapping::CoordResult>>)&&(
                    is_contravariant_vector_index_set_v<
                            ddc::to_type_seq_t<typename Mapping::CoordResult>>),
            "The metric tensor can only be calculated from the Jacobian for coordinate "
            "transformations which map to the Cartesian space.");

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;

public:
    using ContravariantVectorType = DTensor<Dims>;

    using CovariantVectorType = DTensor<vector_index_set_dual_t<Dims>>;

    using CoordArg = PositionCoordinate;

private:
    Mapping m_mapping;

public:
    explicit KOKKOS_FUNCTION MetricTensorEvaluator(Mapping mapping) : m_mapping(mapping) {}

    KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov> operator()(CoordArg const& coord) const
    {
        Tensor J = m_mapping.jacobian_matrix(coord);
        return tensor_mul(index<'j', 'i'>(J), index<'j', 'k'>(J));
    }

    KOKKOS_FUNCTION DTensor<Dims, Dims> inverse(CoordArg const& coord) const
    {
        InverseJacobianMatrix<Mapping, PositionCoordinate> get_inverse_jacobian(m_mapping);
        Tensor inv_J = get_inverse_jacobian(coord);
        return tensor_mul(index<'i', 'j'>(inv_J), index<'k', 'j'>(inv_J));
    }
};
```


