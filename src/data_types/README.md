# Utility Functions

## VectorField

The classes VectorField and VectorFieldSpan provide a way to represent a vector field. The two classes are required as VectorField contains ddc::Chunk objects while VectorFieldSpan contains ddc::ChunkSpan objects.

A VectorField or VectorFieldSpan contains an array of ddc::Chunk or ddc::ChunkSpan objects. When the VectorField is indexed using the `operator()` a ddc::detail::TaggedVector type is returned. This is similiar to a `ddc::Coordinate` and can be accessed in the same way.
Internally the information for the vector fields along each dimension is stored in separate memory chunks. As a result it is not possible to obtain a reference to an element of a VectorField. Instead you must access the internal chunk using `ddcHelper::get`.

In order to facilitate the usage of VectorField the utility functions `ddcHelper::get_domain`, `ddcHelper::deepcopy`, and `ddcHelper::get` are provided. However for the best access it is advised to retrieve the ChunkSpan and use this directly.


## DerivField

The classes DerivField and DerivFieldSpan provide a way to represent a field with its associated derivatives. The two classes are required as DerivField contains ddc::Chunk objects while DerivFieldSpan contains ddc::ChunkSpan objects.

The values and the different combinations of derivatives are each stored in their own chunk as these objects usually have different domains. The chunk itself can be accessed using slicing methods (`operator[]`). The values additionally have a helper method `get_values_span`.
When slicing a DerivField, you can use either domains or elements. These describe the derivatives which should appear in the chunk of interest and any additional information necessary to obtain a ddc::ChunkSpan. If a dimension is not described then it is assumed that the derivative in this direction is not of interest.

Beware: A ddc::ChunkSpan cannot store data defined on a non-contiguous domain (e.g. a DiscreteSubDomain) so when accessing derivatives the position of the derivative must also be included in the slice index.

As for VectorField it is advised to use this object for storage and to interact with the underlying chunks directly. However the utility function `ddcHelper::deepcopy` is nevertheless provided.
