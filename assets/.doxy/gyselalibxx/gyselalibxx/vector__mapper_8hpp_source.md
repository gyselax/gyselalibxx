

# File vector\_mapper.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**vector\_mapper.hpp**](vector__mapper_8hpp.md)

[Go to the documentation of this file](vector__mapper_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

template <class InVectorSpace, class OutVectorSpace, class Mapping, class ExecSpace>
class VectorMapper;

template <class XIn, class YIn, class XOut, class YOut, class Mapping, class ExecSpace>
class VectorMapper<VectorIndexSet<XIn, YIn>, VectorIndexSet<XOut, YOut>, Mapping, ExecSpace>
{
    static_assert(is_accessible_v<ExecSpace, Mapping>);
    static_assert(
            (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordArg>)
            || (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordResult>));
    static_assert(
            (std::is_same_v<Coord<XOut, YOut>, typename Mapping::CoordArg>)
            || (std::is_same_v<Coord<XOut, YOut>, typename Mapping::CoordResult>));

public:
    using memory_space = typename ExecSpace::memory_space;

private:
    Mapping m_mapping;

public:
    explicit VectorMapper(Mapping mapping) : m_mapping(mapping) {}

    template <class IdxRangeType, class LayoutStridedPolicy1, class LayoutStridedPolicy2>
    void operator()(
            ExecSpace exec_space,
            VectorField<
                    double,
                    IdxRangeType,
                    VectorIndexSet<XOut, YOut>,
                    memory_space,
                    LayoutStridedPolicy1> vector_field_output,
            VectorConstField<
                    double,
                    IdxRangeType,
                    VectorIndexSet<XIn, YIn>,
                    memory_space,
                    LayoutStridedPolicy2> vector_field_input)
    {
        using IdxType = typename IdxRangeType::discrete_element_type;

        if constexpr (std::is_same_v<Coord<XIn, YIn>, typename Mapping::CoordArg>) {
            Mapping mapping_proxy = m_mapping;
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(vector_field_input),
                    KOKKOS_LAMBDA(IdxType idx) {
                        Tensor jacobian = mapping_proxy.jacobian_matrix(ddc::coordinate(idx));
                        DVector<XOut, YOut> vector_field_out = tensor_mul(
                                index<'i', 'j'>(jacobian),
                                index<'j'>(vector_field_input));
                        ddcHelper::get<XOut>(vector_field_output)(idx)
                                = ddcHelper::get<XOut>(vector_field_out);
                        ddcHelper::get<YOut>(vector_field_output)(idx)
                                = ddcHelper::get<YOut>(vector_field_out);
                    });
        } else {
            InverseJacobianMatrix<Mapping, ddc::coordinate_of_t<IdxType>> inv_mapping(m_mapping);
            ddc::parallel_for_each(
                    exec_space,
                    get_idx_range(vector_field_input),
                    KOKKOS_LAMBDA(IdxType idx) {
                        Tensor map_J = inv_mapping(ddc::coordinate(idx));

                        DVector<XOut, YOut> vector_out = tensor_mul(
                                index<'i', 'j'>(map_J),
                                index<'j'>(vector_field_input(idx)));
                        ddcHelper::get<XOut>(vector_field_output)(idx)
                                = ddcHelper::get<XOut>(vector_out);
                        ddcHelper::get<YOut>(vector_field_output)(idx)
                                = ddcHelper::get<YOut>(vector_out);
                    });
        }
    }
};


template <
        class ExecSpace,
        class Mapping,
        class ElementType,
        class IdxRangeType,
        class X,
        class Y,
        class LayoutStridedPolicy>
auto create_geometry_mirror_view(
        ExecSpace exec_space,
        VectorField<
                ElementType,
                IdxRangeType,
                VectorIndexSet<X, Y>,
                typename ExecSpace::memory_space,
                LayoutStridedPolicy> vector_field,
        Mapping mapping)
{
    using CoordOut = std::conditional_t<
            std::is_same_v<typename Mapping::CoordArg, Coord<X, Y>>,
            typename Mapping::CoordResult,
            typename Mapping::CoordArg>;
    if constexpr (std::is_same_v<CoordOut, Coord<X, Y>>) {
        return vector_field;
    } else {
        using X_out = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordOut>>;
        using Y_out = ddc::type_seq_element_t<1, ddc::to_type_seq_t<CoordOut>>;
        VectorFieldMem<
                std::remove_const_t<ElementType>,
                IdxRangeType,
                VectorIndexSet<X_out, Y_out>,
                typename ExecSpace::memory_space>
                vector_field_out(get_idx_range(vector_field));
        VectorMapper<VectorIndexSet<X, Y>, VectorIndexSet<X_out, Y_out>, Mapping, ExecSpace>
                vector_mapping(mapping);
        vector_mapping(exec_space, get_field(vector_field_out), get_const_field(vector_field));
        return vector_field_out;
    }
}
```


