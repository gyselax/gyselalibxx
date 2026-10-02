# DDC Quick Reference

This page summarises the [DDC](https://ddc.mdls.fr/) types and functions that are used most often in Gyselalib++.
It is intended as a reference to look things up while writing code.
For an explanation of what each of these types represents and why they are needed, see [DDC in Gyselalib++](DDC_in_gyselalibxx.md).

## DDC and Gyselalib++ names

The names used in DDC are not always intuitive for mathematicians and physicists, so Gyselalib++ provides aliases with more descriptive names.
Gyselalib++ code should use the Gyselalib++ names, however the DDC names will appear in compiler error messages and in the DDC documentation.
The tables below can be used to translate between the two.

### Types

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

### Functions

These functions are defined in `src/utils/ddc_alias_inline_functions.hpp`.
Unlike the DDC methods, they also work with Gyselalib++ types which contain fields (e.g. `VectorField` or `DerivField`).

| Gyselalib++ | DDC | Description |
| --- | --- | --- |
| `get_idx_range(field)` | `field.domain()` | Get the index range on which a field is defined. |
| `get_field(field_mem)` | `field_mem.span_view()` | Get a modifiable `Field` from a `FieldMem` (or a `Field`). |
| `get_const_field(field_mem)` | `field_mem.span_cview()` | Get a `ConstField` from a `FieldMem` (or a `Field`). |
| `get_spline_idx_range(builder)` | `builder.spline_domain()` | Get the index range of the B-splines used by a spline builder. |

### Vocabulary

The same naming differences are found in the vocabulary used in the documentation of the two projects.

| Gyselalib++ | DDC |
| --- | --- |
| Dimension | Continuous dimension |
| Grid | Discrete dimension |
| Index range | Discrete domain |
| Field | Chunk span |

## Common operations

The examples below use a grid `GridX` along a dimension `X`, and the aliases (`IdxX`, `IdxStepX`, `IdxRangeX`, `DFieldMemX`, ...) which are usually defined in the `geometry.hpp` file of each geometry.

### Coordinates and indices

```cpp
IdxX ix(idx_range_x.front());               // The first index of an index range
CoordX x = ddc::coordinate(ix);             // The coordinate of a grid point
double x_value = ddc::get<X>(coord_xy);     // A component of a multi-dimensional coordinate as a double
IdxX ix_next = ix + IdxStepX(1);            // The next point on the grid
IdxStepX n_steps = ix_last - ix_first;      // The number of steps between two indices
double dx = ddc::step<GridX>();             // The distance between two points of a uniform grid

CoordX x_only(coord_xy);                    // The coordinate along one dimension of a multi-dimensional coordinate
IdxX ix_from_xy(ixy);                       // The index along one dimension of a multi-dimensional index
IdxStepX n_steps_x(n_steps_xy);             // The step along one dimension of a multi-dimensional step
IdxRangeX idx_range_x(idx_range_xy);        // The index range along one dimension of a multi-dimensional index range
```

To extract part of a multi-dimensional `Coord`, `Idx`, `IdxStep` or `IdxRange`, use the constructor of the lower-dimensional type as shown above.
This is preferred to the equivalent but more verbose `ddc::select<GridX>(ixy)`.

### Index ranges

```cpp
IdxRangeXY idx_range_xy = get_idx_range(field_xy); // The index range on which a field is defined
IdxRangeX idx_range_x = get_idx_range<GridX>(field_xy); // The index range along one dimension
IdxX ix_first = idx_range_x.front();               // The first index
IdxX ix_last = idx_range_x.back();                 // The last index
std::size_t n_points = idx_range_xy.size();        // The total number of points
IdxStep<GridX, GridY> n_points_per_dim = idx_range_xy.extents(); // The number of points in each dimension
```

An index range covering part of another index range can be created with:

- `idx_range.take_first(n)` : the first `n` indices.
- `idx_range.take_last(n)` : the last `n` indices.
- `idx_range.remove_first(n)` : all but the first `n` indices.
- `idx_range.remove_last(n)` : all but the last `n` indices.
- `idx_range.remove(n_first, n_last)` : all but the first `n_first` and the last `n_last` indices.

### Allocating and accessing data

```cpp
DFieldMemX field_alloc(idx_range_x);                   // Allocate memory (on the GPU by default)
DFieldX field = get_field(field_alloc);                // A modifiable reference to the data
DConstFieldX field_const = get_const_field(field_alloc); // A read-only reference to the data
host_t<DFieldMemX> field_host_alloc(idx_range_x);      // Allocate memory on the CPU
double value = field(ix);                              // Access an element
DFieldY field_y = field_xy[ix];                        // A slice of a multi-dimensional field
DFieldX field_part = field[idx_range_sub];             // A reference to part of a field
```

Elements of a field stored on the GPU can only be accessed inside a GPU kernel (see [Loops](#loops)).
When writing a function, take `Field` or `ConstField` arguments rather than `FieldMem` (see [Parameter passing](../contributing/CODING_STANDARD.md#parameter-passing)).

### Loops

```cpp
// The location of the current function, used to name the parallel loops
std::source_location const location = std::source_location::current();

// Loop over an index range on the GPU
ddc::parallel_for_each(
        location.function_name(),
        Kokkos::DefaultExecutionSpace(),
        idx_range_x,
        KOKKOS_LAMBDA(IdxX ix) { field(ix) = ...; });

// Loop over an index range on the CPU
ddc::host_for_each(idx_range_x, [&](IdxX ix) { field_host(ix) = ...; });

// Compute a reduction (here the maximum value) on the GPU
double max_value = ddc::parallel_transform_reduce(
        location.function_name(),
        Kokkos::DefaultExecutionSpace(),
        idx_range_x,
        0.0,
        ddc::reducer::max<double>(),
        KOKKOS_LAMBDA(IdxX ix) { return field(ix); });
```

Other reducers include `ddc::reducer::sum` and `ddc::reducer::land` (logical and).

Always pass `location.function_name()` as the first argument of `ddc::parallel_for_each` and `ddc::parallel_transform_reduce`. This name identifies the loop in the output of [profiling tools](../how_to/profiling.md).

Always pass an execution space to parallel loops (see [Synchronicity](DDC_in_gyselalibxx.md#synchronicity)).
If the loop does not compile, see [Compilation errors](../troubleshooting/index.md#compilation-errors).

### Copying data

```cpp
ddc::parallel_fill(field, 0.0);                     // Set every element to a value
ddc::parallel_deepcopy(field_dest, field_src);      // Copy data between two fields with the same index range

// Get a copy of GPU data on the CPU (e.g. to save it to a file)
auto field_host = ddc::create_mirror_view_and_copy(field);

// Get a copy of CPU data on the GPU
auto field_device = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), field_host);
```

`ddc::create_mirror_view_and_copy` only allocates new memory if the data is not already accessible from the target memory space.
