

# File curl.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**curl.hpp**](curl_8hpp.md)

[Go to the documentation of this file](curl_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gradient.hpp"
#include "metric_tensor_evaluator.hpp"
#include "static_tensors.hpp"
#include "vector_field.hpp"

template <concepts::Mapping Mapping3D, class MappingCoord = typename Mapping3D::CoordArg>
class Curl
{
    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D m_mapping;
    MetricTensorEvaluator<Mapping3D, MappingCoord> m_metric_tensor;

public:
    explicit Curl(Mapping3D const& mapping) : m_mapping(mapping), m_metric_tensor(mapping) {}

    template <class T, class CovBasisF>
    KOKKOS_INLINE_FUNCTION Tensor<T, get_contravariant_dims_t<CovBasisF>> operator()(
            Tensor<T, CovBasisF, CovBasisSpatial> const& partial_derivatives_f,
            MappingCoord const& coord) const
    {
        static_assert(
                is_covariant_vector_index_set_v<CovBasisF>,
                "Expected the derivative of f expressed in the covariant basis");
        double J = m_mapping.jacobian(coord);
        LeviCivitaTensor<double, BasisSpatial> eps(J);

        // The following code computes: 1/sqrt(g) * eps^{klm} \nabla_{l} F_{m}
        return tensor_mul(index<'k', 'l', 'm'>(eps), index<'m', 'l'>(partial_derivatives_f));
    }
};
```


