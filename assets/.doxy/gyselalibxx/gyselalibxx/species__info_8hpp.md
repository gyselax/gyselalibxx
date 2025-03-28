

# File species\_info.hpp



[**FileList**](files.md) **>** [**speciesinfo**](dir_661be8452a62f1b4720eb6eb57123ae7.md) **>** [**species\_info.hpp**](species__info_8hpp.md)

[Go to the source code of this file](species__info_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Species**](structSpecies.md) <br> |
| class | [**SpeciesInformation**](classSpeciesInformation.md) <br>[_**Species**_](structSpecies.md) _discrete dimension to access constant attributes related to species._ |
| class | [**Impl**](classSpeciesInformation_1_1Impl.md) &lt;class Grid1D, class MemorySpace&gt;<br>[_**Impl**_](classSpeciesInformation_1_1Impl.md) _object storing attributes in_`MemorySpace` _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; ElementType, IdxRangeSp &gt; | [**ConstFieldSp**](#typedef-constfieldsp)  <br> |
| typedef ConstFieldSp&lt; double &gt; | [**DConstFieldSp**](#typedef-dconstfieldsp)  <br> |
| typedef FieldMemSp&lt; double &gt; | [**DFieldMemSp**](#typedef-dfieldmemsp)  <br> |
| typedef FieldSp&lt; double &gt; | [**DFieldSp**](#typedef-dfieldsp)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSp &gt; | [**FieldMemSp**](#typedef-fieldmemsp)  <br> |
| typedef Field&lt; ElementType, IdxRangeSp &gt; | [**FieldSp**](#typedef-fieldsp)  <br> |
| typedef FieldMemSp&lt; int &gt; | [**IFieldMemSp**](#typedef-ifieldmemsp)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md) &gt; | [**IdxRangeSp**](#typedef-idxrangesp)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md) &gt; | [**IdxSp**](#typedef-idxsp)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md) &gt; | [**IdxStepSp**](#typedef-idxstepsp)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_species\_information\_v**](#variable-is_species_information_v)   = `std::is\_same\_v&lt;typename Grid1D::discrete\_dimension\_type, [**SpeciesInformation**](classSpeciesInformation.md)&gt;`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**charge**](#function-charge) (Idx&lt; [**Species**](structSpecies.md) &gt; const isp) <br> |
|  KOKKOS\_INLINE\_FUNCTION Idx&lt; [**Species**](structSpecies.md) &gt; | [**ielec**](#function-ielec) () <br> |
|  KOKKOS\_INLINE\_FUNCTION double | [**mass**](#function-mass) (Idx&lt; [**Species**](structSpecies.md) &gt; const isp) <br> |




























## Public Types Documentation




### typedef ConstFieldSp 

```C++
using ConstFieldSp =  ConstField<ElementType, IdxRangeSp>;
```




<hr>



### typedef DConstFieldSp 

```C++
using DConstFieldSp =  ConstFieldSp<double>;
```




<hr>



### typedef DFieldMemSp 

```C++
using DFieldMemSp =  FieldMemSp<double>;
```




<hr>



### typedef DFieldSp 

```C++
using DFieldSp =  FieldSp<double>;
```




<hr>



### typedef FieldMemSp 

```C++
using FieldMemSp =  FieldMem<ElementType, IdxRangeSp>;
```




<hr>



### typedef FieldSp 

```C++
using FieldSp =  Field<ElementType, IdxRangeSp>;
```




<hr>



### typedef IFieldMemSp 

```C++
using IFieldMemSp =  FieldMemSp<int>;
```




<hr>



### typedef IdxRangeSp 

```C++
using IdxRangeSp =  IdxRange<Species>;
```




<hr>



### typedef IdxSp 

```C++
using IdxSp =  Idx<Species>;
```




<hr>



### typedef IdxStepSp 

```C++
using IdxStepSp =  IdxStep<Species>;
```




<hr>
## Public Attributes Documentation




### variable is\_species\_information\_v 

```C++
constexpr bool is_species_information_v;
```




<hr>
## Public Functions Documentation




### function charge 

```C++
KOKKOS_INLINE_FUNCTION double charge (
    Idx< Species > const isp
) 
```





**Parameters:**


* `isp` a discrete element of either a kinetic or adiabatic species 



**Returns:**

the charge associated to the discrete element 





        

<hr>



### function ielec 

```C++
KOKKOS_INLINE_FUNCTION Idx< Species > ielec () 
```





**Returns:**

the discrete element representing the electron species 





        

<hr>



### function mass 

```C++
KOKKOS_INLINE_FUNCTION double mass (
    Idx< Species > const isp
) 
```





**Parameters:**


* `isp` a discrete element of either a kinetic or adiabatic species 



**Returns:**

the mass associated to the discrete element 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/speciesinfo/species_info.hpp`

