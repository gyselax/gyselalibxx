

# File tensor.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor.hpp**](tensor_8hpp.md)

[Go to the source code of this file](tensor_8hpp_source.md)



* `#include <array>`
* `#include <cassert>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "tensor_index_tools.hpp"`
* `#include "vector_index_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**Tensor**](classTensor.md) &lt;class ElementType, ValidIndexSet&gt;<br>_A class representing a_ [_**Tensor**_](classTensor.md) _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Tensor**](classTensor.md)&lt; double, ValidIndexSet... &gt; | [**DTensor**](#typedef-dtensor)  <br>_A helper type alias to get a tensor containing doubles._  |
| typedef [**Vector**](classTensor.md)&lt; double, Dims... &gt; | [**DVector**](#typedef-dvector)  <br>_A helper type alias to get a 1D tensor (a vector) of doubles._  |
| typedef [**Tensor**](classTensor.md)&lt; ElementType, VectorIndexSet&lt; Dims... &gt; &gt; | [**Vector**](#typedef-vector)  <br>_A helper type alias to get a 1D tensor (a vector)._  |
| typedef typename detail::ToTensor&lt; ElementType, TypeSeqValidIndexSet &gt;::type | [**to\_tensor\_t**](#typedef-to_tensor_t)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, ValidIndexSet... &gt; | [**operator\***](#function-operator) (OElementType val, [**Tensor**](classTensor.md)&lt; ElementType, ValidIndexSet... &gt; const & tensor) <br>_An operator to multiply all a tensor by a value._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; | [**operator+**](#function-operator_1) (Coord&lt; Dims... &gt; const & coord, [**DVector**](classTensor.md)&lt; Dims... &gt; const & tensor) <br> |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; & | [**operator+=**](#function-operator_2) (Coord&lt; Dims... &gt; & coord, [**DVector**](classTensor.md)&lt; Dims... &gt; const & tensor) <br> |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; | [**operator-**](#function-operator-) (Coord&lt; Dims... &gt; const & coord, [**DVector**](classTensor.md)&lt; Dims... &gt; const & tensor) <br> |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; & | [**operator-=**](#function-operator-_1) (Coord&lt; Dims... &gt; & coord, [**DVector**](classTensor.md)&lt; Dims... &gt; const & tensor) <br> |




























## Public Types Documentation




### typedef DTensor 

_A helper type alias to get a tensor containing doubles._ 
```C++
using DTensor =  Tensor<double, ValidIndexSet...>;
```





**Template parameters:**


* `ValidIndexSet` The indices that can be used along each dimension of the tensor. 




        

<hr>



### typedef DVector 

_A helper type alias to get a 1D tensor (a vector) of doubles._ 
```C++
using DVector =  Vector<double, Dims...>;
```





**Template parameters:**


* `Dims` The dimensions that can be used to index the vector. 




        

<hr>



### typedef Vector 

_A helper type alias to get a 1D tensor (a vector)._ 
```C++
using Vector =  Tensor<ElementType, VectorIndexSet<Dims...> >;
```





**Template parameters:**


* `ElementType` The type of the elements of the tensor (usually double/complex). 
* `Dims` The dimensions that can be used to index the vector. 




        

<hr>



### typedef to\_tensor\_t 

```C++
using to_tensor_t =  typename detail::ToTensor<ElementType, TypeSeqValidIndexSet>::type;
```




<hr>
## Public Functions Documentation




### function operator\* 

_An operator to multiply all a tensor by a value._ 
```C++
template<class ElementType, class OElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, ValidIndexSet... > operator* (
    OElementType val,
    Tensor < ElementType, ValidIndexSet... > const & tensor
) 
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 
* `tensor` The tensor being multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>



### function operator+ 

```C++
template<class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > operator+ (
    Coord< Dims... > const & coord,
    DVector < Dims... > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate to which the tensor is added. 
* `tensor` The tensor to be added to the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator+= 

```C++
template<class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > & operator+= (
    Coord< Dims... > & coord,
    DVector < Dims... > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate to which the tensor is added. 
* `tensor` The tensor to be added to the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator- 

```C++
template<class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > operator- (
    Coord< Dims... > const & coord,
    DVector < Dims... > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate from which the tensor is subtracted. 
* `tensor` The tensor to be subtracted from the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator-= 

```C++
template<class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > & operator-= (
    Coord< Dims... > & coord,
    DVector < Dims... > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate from which the tensor is subtracted. 
* `tensor` The tensor to be subtracted from the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor.hpp`

