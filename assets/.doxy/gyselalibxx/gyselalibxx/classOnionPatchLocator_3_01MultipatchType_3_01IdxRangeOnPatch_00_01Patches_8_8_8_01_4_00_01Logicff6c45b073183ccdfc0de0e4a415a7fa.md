

# Class OnionPatchLocator&lt; MultipatchType&lt; IdxRangeOnPatch, Patches... &gt;, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace &gt;

**template &lt;class... Patches, class LogicalToPhysicalMapping, class PhysicalToLogicalMapping, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**OnionPatchLocator&lt; MultipatchType&lt; IdxRangeOnPatch, Patches... &gt;, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace &gt;**](classOnionPatchLocator_3_01MultipatchType_3_01IdxRangeOnPatch_00_01Patches_8_8_8_01_4_00_01Logicff6c45b073183ccdfc0de0e4a415a7fa.md)



[_**Patch**_](structPatch.md) _locator specialised for "onion" geometry._[More...](#detailed-description)

* `#include <onion_patch_locator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**MultipatchType**](classMultipatchType.md)&lt; IdxRangeOnPatch, Patches... &gt; | [**MultipatchIdxRanges**](#typedef-multipatchidxranges)  <br>[_**MultipatchType**_](classMultipatchType.md) _storing the index range._ |
| typedef ddc::detail::TypeSeq&lt; Patches... &gt; | [**PatchOrdering**](#typedef-patchordering)  <br>_Sequence ddc::detail::TypeSeq of patch tags._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |
| typedef LogicalToPhysicalMapping | [**get\_mapping\_on\_logical\_dim\_t**](#typedef-get_mapping_on_logical_dim_t)  <br>_Get the type of the mapping from given logical continuous dimensions._  |
| typedef LogicalToPhysicalMapping | [**get\_mapping\_on\_patch\_t**](#typedef-get_mapping_on_patch_t)  <br>_Get the type of the mapping on the given_ [_**Patch**_](structPatch.md) _._ |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr int | [**outside\_rmax\_domain**](#variable-outside_rmax_domain)   = `-1`<br>_Default value to define outside domain (not a patch) for radius bigger than the maximum radius._  |
|  constexpr int | [**outside\_rmin\_domain**](#variable-outside_rmin_domain)   = `-2`<br>_Default value to define outside domain (not a patch) for radius smaller than the minimum radius._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**OnionPatchLocator**](#function-onionpatchlocator) ([**MultipatchIdxRanges**](classOnionPatchLocator_3_01MultipatchType_3_01IdxRangeOnPatch_00_01Patches_8_8_8_01_4_00_01Logicff6c45b073183ccdfc0de0e4a415a7fa.md#typedef-multipatchidxranges) const & all\_idx\_ranges, LogicalToPhysicalMapping const & to\_physical\_mapping, PhysicalToLogicalMapping const & to\_logical\_mapping) <br>_Instantiante the operator with_ [_**MultipatchType**_](classMultipatchType.md) _of index ranges and a mapping on all the patches._ |
|  KOKKOS\_FUNCTION LogicalToPhysicalMapping | [**get\_mapping\_on\_logical\_dim**](#function-get_mapping_on_logical_dim) () const<br>_Get the mapping from given logical continuous dimensions. The function can run on device and host._  |
|  KOKKOS\_FUNCTION LogicalToPhysicalMapping | [**get\_mapping\_on\_patch**](#function-get_mapping_on_patch) () const<br>_Get the mapping on the given_ [_**Patch**_](structPatch.md) _. The function can run on device and host._ |
|  KOKKOS\_INLINE\_FUNCTION int | [**operator()**](#function-operator) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const coord) const<br>_Get the patch where the given physical coordinate is._  |
|   | [**~OnionPatchLocator**](#function-onionpatchlocator) () = default<br> |




























## Detailed Description


See OnionPatchLocatorImplementation




**Template parameters:**


* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 
* `Patches` [**Patch**](structPatch.md) types. Their order is important. 
* `LogicalToPhysicalMapping` A mapping type for all the patches. 




    
## Public Types Documentation




### typedef MultipatchIdxRanges 

[_**MultipatchType**_](classMultipatchType.md) _storing the index range._
```C++
using OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::MultipatchIdxRanges =  MultipatchType<IdxRangeOnPatch, Patches...>;
```




<hr>



### typedef PatchOrdering 

_Sequence ddc::detail::TypeSeq of patch tags._ 
```C++
using OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::PatchOrdering =  ddc::detail::TypeSeq<Patches...>;
```




<hr>



### typedef exec\_space 

_The space (CPU/GPU) where the calculations are carried out._ 
```C++
using OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::exec_space =  ExecSpace;
```




<hr>



### typedef get\_mapping\_on\_logical\_dim\_t 

_Get the type of the mapping from given logical continuous dimensions._ 
```C++
using OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::get_mapping_on_logical_dim_t =  LogicalToPhysicalMapping;
```




<hr>



### typedef get\_mapping\_on\_patch\_t 

_Get the type of the mapping on the given_ [_**Patch**_](structPatch.md) _._
```C++
using OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::get_mapping_on_patch_t =  LogicalToPhysicalMapping;
```




<hr>
## Public Static Attributes Documentation




### variable outside\_rmax\_domain 

_Default value to define outside domain (not a patch) for radius bigger than the maximum radius._ 
```C++
constexpr int OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::outside_rmax_domain;
```




<hr>



### variable outside\_rmin\_domain 

_Default value to define outside domain (not a patch) for radius smaller than the minimum radius._ 
```C++
constexpr int OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::outside_rmin_domain;
```




<hr>
## Public Functions Documentation




### function OnionPatchLocator 

_Instantiante the operator with_ [_**MultipatchType**_](classMultipatchType.md) _of index ranges and a mapping on all the patches._
```C++
inline OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::OnionPatchLocator (
    MultipatchIdxRanges const & all_idx_ranges,
    LogicalToPhysicalMapping const & to_physical_mapping,
    PhysicalToLogicalMapping const & to_logical_mapping
) 
```



The order of the elements in the tuple or the [**MultipatchType**](classMultipatchType.md) doesn't matter.




**Parameters:**


* `all_idx_ranges` A [**MultipatchType**](classMultipatchType.md) of index ranges defined on the logical domain of each patch. 
* `to_physical_mapping` Mapping from the logical domains of every patch to the global physical domain. 
* `to_logical_mapping` Mapping from the global physical domain to the logical domain of every patch. 




        

<hr>



### function get\_mapping\_on\_logical\_dim 

_Get the mapping from given logical continuous dimensions. The function can run on device and host._ 
```C++
template<class Dim1, class Dim2>
inline KOKKOS_FUNCTION LogicalToPhysicalMapping OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::get_mapping_on_logical_dim () const
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 



**Returns:**

The mapping on the given [**Patch**](structPatch.md). 





        

<hr>



### function get\_mapping\_on\_patch 

_Get the mapping on the given_ [_**Patch**_](structPatch.md) _. The function can run on device and host._
```C++
template<class Patch>
inline KOKKOS_FUNCTION LogicalToPhysicalMapping OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::get_mapping_on_patch () const
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 



**Returns:**

The mapping on the given [**Patch**](structPatch.md). 





        

<hr>



### function operator() 

_Get the patch where the given physical coordinate is._ 
```C++
inline KOKKOS_INLINE_FUNCTION int OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::operator() (
    Coord< X , Y > const coord
) const
```



We use a dichotomy method to find the patch the physical coordinate is on.


Each logical grid of the patches are defined on the same dimensions. Knowing that, we can compare the logical coordinates between the patches in the dichotomy.




**Parameters:**


* `coord` [in] The given physical coordinate. 



**Returns:**

[int] The patch index where the physical coordinate. If the coordinate is outside of the domain, it returns a negative value. 





        

<hr>



### function ~OnionPatchLocator 

```C++
OnionPatchLocator< MultipatchType< IdxRangeOnPatch, Patches... >, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace >::~OnionPatchLocator () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/onion_patch_locator.hpp`

