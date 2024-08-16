# Data Storage Types

## VectorField

The classes VectorFieldMem and VectorField provide a way to represent a vector field. The two classes are required as VectorFieldMem contains FieldMem objects while VectorField contains Field objects.

A VectorFieldMem or VectorField contains an array of FieldMem or Field objects. When the VectorField is indexed using the `operator()` a ddc::detail::TaggedVector type is returned. This is similiar to a `ddc::Coordinate` and can be accessed in the same way.
Internally the information for the vector fields along each dimension is stored in separate memory chunks. As a result it is not possible to obtain a reference to an element of a VectorField. Instead you must access the internal chunk using `ddcHelper::get`.

In order to facilitate the usage of VectorField the utility functions `get_idx_range`, `ddcHelper::deepcopy`, and `ddcHelper::get` are provided. However for the best access it is advised to retrieve the Field and use this directly.


## DerivField

The classes DerivFieldMem and DerivField provide a way to represent a field with its associated derivatives. The two classes are required as DerivFieldMem contains FieldMem objects while DerivField contains Field objects.

The values and the different combinations of derivatives are each stored in their own field as these objects usually have different index ranges. The field itself can be accessed using slicing methods (`operator[]`). The values additionally have a helper method `get_values_field`.
When slicing a DerivField, you can use either idx ranges or individual indices. These describe the derivatives which should appear in the field of interest and any additional information necessary to obtain a Field. If a dimension is not described then it is assumed that the derivative in this direction is not of interest.

Beware: A DDC Field cannot store data defined on a non-contiguous index range (e.g. an IdxRangeSlice) so when accessing derivatives the position of the derivative must also be included in the slice index.

As for VectorField it is advised to use this object for storage and to interact with the underlying fields directly. However the utility function `ddcHelper::deepcopy` is nevertheless provided.
