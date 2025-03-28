

# File ddc\_helper.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_helper.hpp**](ddc__helper_8hpp.md)

[Go to the source code of this file](ddc__helper_8hpp_source.md)



* `#include <cassert>`
* `#include <cmath>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |




## Public Types

| Type | Name |
| ---: | :--- |
| typedef on\_memory\_space\_t&lt; Kokkos::DefaultExecutionSpace::memory\_space, C &gt; | [**device\_t**](#typedef-device_t)  <br>_Alias template helper returning the "device" version of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _._ |
| typedef on\_memory\_space\_t&lt; Kokkos::HostSpace, C &gt; | [**host\_t**](#typedef-host_t)  <br>_Alias template helper returning the "host" version of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _._ |
| typedef typename detail::OnMemorySpace&lt; MemorySpace, C &gt;::type | [**on\_memory\_space\_t**](#typedef-on_memory_space_t)  <br>_Alias template helper returning the type of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _on a MemorySpace._ |
















































## Public Types Documentation




### typedef device\_t 

_Alias template helper returning the "device" version of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _._
```C++
using device_t =  on_memory_space_t<Kokkos::DefaultExecutionSpace::memory_space, C>;
```




<hr>



### typedef host\_t 

_Alias template helper returning the "host" version of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _._
```C++
using host_t =  on_memory_space_t<Kokkos::HostSpace, C>;
```




<hr>



### typedef on\_memory\_space\_t 

_Alias template helper returning the type of a_ `ddc::Chunk` _, a_`ddc::ChunkSpan` _, a_`VectorFieldMem` _or a_`VectorField` _on a MemorySpace._
```C++
using on_memory_space_t =  typename detail::OnMemorySpace<MemorySpace, C>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/ddc_helper.hpp`

