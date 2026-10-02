# DDC and Gyselalib++ Names

The names used in [DDC](https://ddc.mdls.fr/) are not always intuitive for mathematicians and physicists, so Gyselalib++ provides aliases with more descriptive names.
Gyselalib++ code should use the Gyselalib++ names, however the DDC names will appear in compiler error messages and in the DDC documentation.
This page can be used to translate between the two.
For an explanation of what each of these types represents, see [DDC in Gyselalib++](DDC_in_gyselalibxx.md).

## Types

These aliases are defined in `src/utils/ddc_aliases.hpp`.

| Gyselalib++ | DDC | Description |
| --- | --- | --- |
| `Coord<Dims...>` | `ddc::Coordinate<Dims...>` | A [coordinate](DDC_in_gyselalibxx.md#coordinates) in continuous space. |
| `Idx<Grids...>` | `ddc::DiscreteElement<Grids...>` | An [index](DDC_in_gyselalibxx.md#index) identifying a point on a grid. |
| `IdxStep<Grids...>` | `ddc::DiscreteVector<Grids...>` | The [number of steps](DDC_in_gyselalibxx.md#index-step) between two indices. |
| `IdxRange<Grids...>` | `ddc::DiscreteDomain<Grids...>` | A contiguous [range of indices](DDC_in_gyselalibxx.md#index-range) on which data is defined. |
| `IdxRangeSlice<Grids...>` | `ddc::StridedDiscreteDomain<Grids...>` | A range of indices with a constant stride between them. |
| `FieldMem<T, IdxRange, MemSpace>` | `ddc::Chunk<T, IdxRange, ddc::KokkosAllocator<T, MemSpace>>` | An object which [allocates the memory](DDC_in_gyselalibxx.md#data-storage) for a field. |
| `Field<T, IdxRange, MemSpace, Layout>` | `ddc::ChunkSpan<T, IdxRange, Layout, MemSpace>` | A modifiable [reference to the data](DDC_in_gyselalibxx.md#data-storage) of a field. Copying it does not allocate memory. |
| `ConstField<T, IdxRange, MemSpace, Layout>` | `ddc::ChunkView<T, IdxRange, Layout, MemSpace>` | A read-only reference to the data of a field. |
| `DFieldMem`, `DField`, `DConstField` | | The same as `FieldMem`, `Field` and `ConstField` with `T = double`. |
| `UniformGridBase<Dim>` | `ddc::UniformPointSampling<Dim>` | The base class of a [grid](DDC_in_gyselalibxx.md#grid) whose points are uniformly spaced. |
| `NonUniformGridBase<Dim>` | `ddc::NonUniformPointSampling<Dim>` | The base class of a grid whose points are not uniformly spaced. |

## Functions

These functions are defined in `src/utils/ddc_alias_inline_functions.hpp`.
Unlike the DDC methods, they also work with Gyselalib++ types which contain fields (e.g. `VectorField` or `DerivField`).

| Gyselalib++ | DDC | Description |
| --- | --- | --- |
| `get_idx_range(field)` | `field.domain()` | Get the index range on which a field is defined. |
| `get_field(field_mem)` | `field_mem.span_view()` | Get a modifiable `Field` from a `FieldMem` (or a `Field`). |
| `get_const_field(field_mem)` | `field_mem.span_cview()` | Get a `ConstField` from a `FieldMem` (or a `Field`). |
| `get_spline_idx_range(builder)` | `builder.spline_domain()` | Get the index range of the B-splines used by a spline builder. |

## Vocabulary

The same naming differences are found in the vocabulary used in the documentation of the two projects.

| Gyselalib++ | DDC |
| --- | --- |
| Dimension | Continuous dimension |
| Grid | Discrete dimension |
| Index range | Discrete domain |
| Field | Chunk span |
