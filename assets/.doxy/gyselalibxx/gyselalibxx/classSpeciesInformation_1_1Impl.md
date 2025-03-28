

# Class SpeciesInformation::Impl

**template &lt;class Grid1D, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**SpeciesInformation**](classSpeciesInformation.md) **>** [**Impl**](classSpeciesInformation_1_1Impl.md)



[_**Impl**_](classSpeciesInformation_1_1Impl.md) _object storing attributes in_`MemorySpace` _._

* `#include <species_info.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**SpeciesInformation**](classSpeciesInformation.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_alias of the discrete dimension_  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Impl**](#function-impl-23) ([**Impl**](classSpeciesInformation_1_1Impl.md)&lt; Grid1D, OMemorySpace &gt; const & impl) <br>_Conversion constructor between different memory spaces._  |
|   | [**Impl**](#function-impl-33) (DFieldMem&lt; index\_range\_type, MemorySpace &gt; charge, DFieldMem&lt; index\_range\_type, MemorySpace &gt; mass) <br>_Main constructor taking all attributes._  |
|  KOKKOS\_FUNCTION double | [**charge**](#function-charge) (discrete\_element\_type const isp) const<br> |
|  auto | [**charges**](#function-charges) () const<br> |
|  KOKKOS\_FUNCTION discrete\_element\_type | [**ielec**](#function-ielec) () const<br> |
|  KOKKOS\_FUNCTION double | [**mass**](#function-mass) (discrete\_element\_type const isp) const<br> |
|  auto | [**masses**](#function-masses) () const<br> |




























## Public Types Documentation




### typedef discrete\_dimension\_type 

_alias of the discrete dimension_ 
```C++
using SpeciesInformation::Impl< Grid1D, MemorySpace >::discrete_dimension_type =  SpeciesInformation;
```




<hr>
## Public Functions Documentation




### function Impl [2/3]

_Conversion constructor between different memory spaces._ 
```C++
template<class OMemorySpace>
inline explicit SpeciesInformation::Impl::Impl (
    Impl < Grid1D, OMemorySpace > const & impl
) 
```





**Parameters:**


* `impl` object from `OMemorySpace` that will be used to initialise this object on `MemorySpace` 




        

<hr>



### function Impl [3/3]

_Main constructor taking all attributes._ 
```C++
inline SpeciesInformation::Impl::Impl (
    DFieldMem< index_range_type, MemorySpace > charge,
    DFieldMem< index_range_type, MemorySpace > mass
) 
```





**Parameters:**


* `charge` array storing both kinetic and adiabatic charges 
* `mass` array storing both kinetic and adiabatic masses 




        

<hr>



### function charge 

```C++
inline KOKKOS_FUNCTION double SpeciesInformation::Impl::charge (
    discrete_element_type const isp
) const
```





**Parameters:**


* `isp` a discrete element of either a kinetic or adiabatic species 



**Returns:**

the charge associated to the discrete element 





        

<hr>



### function charges 

```C++
inline auto SpeciesInformation::Impl::charges () const
```





**Returns:**

kinetic and adiabatic charges array 





        

<hr>



### function ielec 

```C++
inline KOKKOS_FUNCTION discrete_element_type SpeciesInformation::Impl::ielec () const
```





**Returns:**

the discrete element representing the electron species 





        

<hr>



### function mass 

```C++
inline KOKKOS_FUNCTION double SpeciesInformation::Impl::mass (
    discrete_element_type const isp
) const
```





**Parameters:**


* `isp` a discrete element of either a kinetic or adiabatic species 



**Returns:**

the mass associated to the discrete element 





        

<hr>



### function masses 

```C++
inline auto SpeciesInformation::Impl::masses () const
```





**Returns:**

kinetic and adiabatic masses array 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/speciesinfo/species_info.hpp`

