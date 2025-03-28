

# File polar\_spline.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**polar\_spline.hpp**](polar__spline_8hpp.md)

[Go to the source code of this file](polar__spline_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_helper.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**ConstPolarSpline**](structConstPolarSpline.md) &lt;class PolarBSplinesType, class MemSpace&gt;<br>_A structure containing the two ConstFields necessary to define a constant reference to a spline on a set of polar basis splines._  |
| struct | [**PolarSpline**](structPolarSpline.md) &lt;class PolarBSplinesType, class MemSpace&gt;<br>_A structure containing the two Fields necessary to define a reference to a spline on a set of polar basis splines._  |
| struct | [**PolarSplineMem**](structPolarSplineMem.md) &lt;class PolarBSplinesType, class MemSpace&gt;<br>_A structure containing the two FieldMems necessary to define a spline on a set of polar basis splines._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_data\_access\_methods&lt; ConstPolarSpline&lt; PolarBSplines, MemorySpace &gt; &gt;**](#variable-enable_data_access_methods-constpolarspline-polarbsplines-memoryspace)   = `true`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; PolarSpline&lt; PolarBSplines, MemorySpace &gt; &gt;**](#variable-enable_data_access_methods-polarspline-polarbsplines-memoryspace)   = `true`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; PolarSplineMem&lt; PolarBSplines, MemorySpace &gt; &gt;**](#variable-enable_data_access_methods-polarsplinemem-polarbsplines-memoryspace)   = `true`<br> |
|  constexpr bool | [**enable\_mem\_type&lt; PolarSplineMem&lt; PolarBSplines, MemorySpace &gt; &gt;**](#variable-enable_mem_type-polarsplinemem-polarbsplines-memoryspace)   = `true`<br> |
|  constexpr bool | [**is\_polar\_spline\_v**](#variable-is_polar_spline_v)   = `false`<br> |
|  constexpr bool | [**is\_polar\_spline\_v&lt; ConstPolarSpline&lt; PolarBSplinesType, MemSpace &gt; &gt;**](#variable-is_polar_spline_v-constpolarspline-polarbsplinestype-memspace)   = `true`<br> |
|  constexpr bool | [**is\_polar\_spline\_v&lt; PolarSpline&lt; PolarBSplinesType, MemSpace &gt; &gt;**](#variable-is_polar_spline_v-polarspline-polarbsplinestype-memspace)   = `true`<br> |
|  constexpr bool | [**is\_polar\_spline\_v&lt; PolarSplineMem&lt; PolarBSplinesType, MemSpace &gt; &gt;**](#variable-is_polar_spline_v-polarsplinemem-polarbsplinestype-memspace)   = `true`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, typename ExecSpace::memory\_space &gt; | [**create\_mirror**](#function-create_mirror) (ExecSpace const & exec\_space, [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create a_ [_**PolarSplineMem**_](structPolarSplineMem.md) _instance, which has the same attributes as src, except the memory space which is related to ExecSpace._ |
|  [**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, Kokkos::HostSpace &gt; | [**create\_mirror**](#function-create_mirror) ([**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create a_ [_**PolarSplineMem**_](structPolarSplineMem.md) _instance, which has the same attributes as src, This function allocates memory on the host._ |
|  [**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, typename ExecSpace::memory\_space &gt; | [**create\_mirror\_and\_copy**](#function-create_mirror_and_copy) (ExecSpace const & exec\_space, [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create copies of_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class instances. It first creates an instance on a memory space accessible from the specified execution space, then it copies the data to the new instance._ |
|  [**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, Kokkos::HostSpace &gt; | [**create\_mirror\_and\_copy**](#function-create_mirror_and_copy) ([**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create copies for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class.It first creates an instance of the class on host, then it operates a copy of the data._ |
|  auto | [**create\_mirror\_view**](#function-create_mirror_view) (ExecSpace const & exec\_space, [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create a mirror view instance for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class on the specified execution space._ |
|  auto | [**create\_mirror\_view**](#function-create_mirror_view) ([**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create a host allocated mirror view for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class._ |
|  auto | [**create\_mirror\_view\_and\_copy**](#function-create_mirror_view_and_copy) (ExecSpace const & exec\_space, [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create copies for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class. If src is accessible from exec\_space, src is returned, else, it first creates a mirror view of the class, then it copies the data to the new instance._ |
|  auto | [**create\_mirror\_view\_and\_copy**](#function-create_mirror_view_and_copy) ([**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & src) <br>_A function to create host allocated view for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class._ |




























## Public Attributes Documentation




### variable enable\_data\_access\_methods&lt; ConstPolarSpline&lt; PolarBSplines, MemorySpace &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< ConstPolarSpline< PolarBSplines, MemorySpace > >;
```




<hr>



### variable enable\_data\_access\_methods&lt; PolarSpline&lt; PolarBSplines, MemorySpace &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< PolarSpline< PolarBSplines, MemorySpace > >;
```




<hr>



### variable enable\_data\_access\_methods&lt; PolarSplineMem&lt; PolarBSplines, MemorySpace &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< PolarSplineMem< PolarBSplines, MemorySpace > >;
```




<hr>



### variable enable\_mem\_type&lt; PolarSplineMem&lt; PolarBSplines, MemorySpace &gt; &gt; 

```C++
constexpr bool enable_mem_type< PolarSplineMem< PolarBSplines, MemorySpace > >;
```




<hr>



### variable is\_polar\_spline\_v 

```C++
constexpr bool is_polar_spline_v;
```




<hr>



### variable is\_polar\_spline\_v&lt; ConstPolarSpline&lt; PolarBSplinesType, MemSpace &gt; &gt; 

```C++
constexpr bool is_polar_spline_v< ConstPolarSpline< PolarBSplinesType, MemSpace > >;
```




<hr>



### variable is\_polar\_spline\_v&lt; PolarSpline&lt; PolarBSplinesType, MemSpace &gt; &gt; 

```C++
constexpr bool is_polar_spline_v< PolarSpline< PolarBSplinesType, MemSpace > >;
```




<hr>



### variable is\_polar\_spline\_v&lt; PolarSplineMem&lt; PolarBSplinesType, MemSpace &gt; &gt; 

```C++
constexpr bool is_polar_spline_v< PolarSplineMem< PolarBSplinesType, MemSpace > >;
```




<hr>
## Public Functions Documentation




### function create\_mirror 

_A function to create a_ [_**PolarSplineMem**_](structPolarSplineMem.md) _instance, which has the same attributes as src, except the memory space which is related to ExecSpace._
```C++
template<class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSplineMem < PolarBSplinesType, typename ExecSpace::memory_space > create_mirror (
    ExecSpace const & exec_space,
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `exec_space` Execution space from which the result must be accessible. 
* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror 

_A function to create a_ [_**PolarSplineMem**_](structPolarSplineMem.md) _instance, which has the same attributes as src, This function allocates memory on the host._
```C++
template<class PolarBSplinesType, class MemSpace>
PolarSplineMem < PolarBSplinesType, Kokkos::HostSpace > create_mirror (
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_and\_copy 

_A function to create copies of_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class instances. It first creates an instance on a memory space accessible from the specified execution space, then it copies the data to the new instance._
```C++
template<class ExecSpace, class PolarBSplinesType, class MemSpace>
PolarSplineMem < PolarBSplinesType, typename ExecSpace::memory_space > create_mirror_and_copy (
    ExecSpace const & exec_space,
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `exec_space` Execution space for allocation. 
* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_and\_copy 

_A function to create copies for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class.It first creates an instance of the class on host, then it operates a copy of the data._
```C++
template<class PolarBSplinesType, class MemSpace>
PolarSplineMem < PolarBSplinesType, Kokkos::HostSpace > create_mirror_and_copy (
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_view 

_A function to create a mirror view instance for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class on the specified execution space._
```C++
template<class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view (
    ExecSpace const & exec_space,
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `exec_space` Execution space from which the result must be accessible. 
* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_view 

_A function to create a host allocated mirror view for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class._
```C++
template<class PolarBSplinesType, class MemSpace>
auto create_mirror_view (
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_view\_and\_copy 

_A function to create copies for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class. If src is accessible from exec\_space, src is returned, else, it first creates a mirror view of the class, then it copies the data to the new instance._
```C++
template<class ExecSpace, class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy (
    ExecSpace const & exec_space,
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**Parameters:**


* `exec_space` Execution space for allocation. 
* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>



### function create\_mirror\_view\_and\_copy 

_A function to create host allocated view for_ [_**PolarSplineMem**_](structPolarSplineMem.md) _class._
```C++
template<class PolarBSplinesType, class MemSpace>
auto create_mirror_view_and_copy (
    PolarSpline < PolarBSplinesType, MemSpace > const & src
) 
```





**See also:** create\_mirror\_view\_and\_copy.


**Parameters:**


* `src` A reference to [**PolarSplineMem**](structPolarSplineMem.md). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_spline.hpp`

