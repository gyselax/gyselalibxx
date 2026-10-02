# Core Concepts

The pages in this section explain the ideas that are used throughout Gyselalib++.
Some of them are essential reading for new users, but they are also intended as a reference to return to whenever you are unsure how something works.

## [DDC in Gyselalib++](DDC_in_gyselalibxx.md)

Read this when you need to understand how data is indexed and stored. It describes:

- [Coordinates](DDC_in_gyselalibxx.md#coordinates) (`Coord`) : positions in continuous space.
- [Grids](DDC_in_gyselalibxx.md#grid) : the discretisation of a continuous dimension.
- [Indices](DDC_in_gyselalibxx.md#index) (`Idx`) and [index steps](DDC_in_gyselalibxx.md#index-step) (`IdxStep`) : positions on a grid and the distance between them.
- [Index ranges](DDC_in_gyselalibxx.md#index-range) (`IdxRange`) : the sets of indices on which data is defined.
- [Data storage](DDC_in_gyselalibxx.md#data-storage) (`FieldMem`, `Field`, `ConstField`) : how memory is allocated and accessed.
- [Pitfalls](DDC_in_gyselalibxx.md#pitfalls) : common mistakes, such as forgetting to synchronise CPU and GPU operations.

## [DDC and Gyselalib++ Names](ddc_names.md)

Read this when you need to translate between the names used in DDC (which appear in compiler error messages and in the DDC documentation) and the names used in Gyselalib++, e.g. `ddc::DiscreteElement` and `Idx`.

## [Mathematical and Physical Conventions](mathematical_and_physical_conventions.md)

Read this when you need to know which mathematical conventions are used in the code. It covers curvilinear coordinates, contravariant and covariant bases, the metric tensor, the Jacobian and the definition of differential operators such as the gradient and the divergence.

## [Covariant and Contravariant Tensors](coding_covariant_and_contravariant_tensors.md)

Read this when you are writing code which uses vectors, tensors or vector fields. It explains how the distinction between covariant and contravariant components is encoded in the type system, how tensors are multiplied, and how to change between vector spaces.
It assumes that you are familiar with the [Mathematical and Physical Conventions](mathematical_and_physical_conventions.md).
