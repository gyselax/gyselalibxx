# Data Types for Multipatch Geometry

This directory contains classes which are useful for handling objects and
types defined on a multipatch domain.


## MultipatchType

The class `MultipatchType` stores different objects of a type which is templated with 
patches (see `Patch`). The `MultipatchType::get` method can be used to 
retrieve an object defined on a specified patch.
For example we can use fields on different patches which would be the type
`DField<Patch::IdxRange12>`.

So after defining 
```
template<class Patch>
using DFieldOnPatch = DField<Patch::IdxRange12>;
```
we could then have three fields `field1`, `field2` and `field3` on 
patches 1,2 and 3 repectively. The `MultipatchType` object would then 
be initialized as
```
MultipatchType<DFieldOnPatch, Patch1, Patch2, Patch3> multipatch_field(field1, field2, field3);
```
and the field on patch 3 can be retrieved via
```
DField<Patch3::IdxRange12> field3_from_multipatch = multipatch_field.get<Patch3>();
```


### Types
In `types.hpp` file different aliases are defined to use the MultipatchType class: 

* DFieldOnPatch: for storing Field defined on the logical 2D grid of patches.
* IdxRangeOnPatch: for storing 2D IdxRange defined on the logical 2D grid of patches.