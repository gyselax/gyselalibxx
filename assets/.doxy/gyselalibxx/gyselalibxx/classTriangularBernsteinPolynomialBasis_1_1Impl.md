

# Class TriangularBernsteinPolynomialBasis::Impl

**template &lt;class DDim, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md) **>** [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md)



[More...](#detailed-description)

* `#include <bernstein.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_The tag which identifies the basis._  |
| typedef IdxRange&lt; DDim &gt; | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The type of the index range of the basis._  |
| typedef Idx&lt; DDim &gt; | [**discrete\_element\_type**](#typedef-discrete_element_type)  <br>_The type of an index of an element of the basis._  |
| typedef IdxStep&lt; DDim &gt; | [**discrete\_vector\_type**](#typedef-discrete_vector_type)  <br>_The type of an index step from one element of the basis to another._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Impl**](#function-impl-25) ([**CartesianToBarycentric**](classCartesianToBarycentric.md)&lt; [**X**](structX.md), [**Y**](structY.md), Corner1Tag, Corner2Tag, Corner3Tag &gt; const & coord\_changer) <br>_Construct the basis from the barycentric coordinate mapping._  |
|   | [**Impl**](#function-impl-35) ([**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md)&lt; DDim, OriginMemorySpace &gt; const & impl) <br>_Construct the basis by copy. This constructor is used to create the class on a different memory space._  |
|   | [**Impl**](#function-impl-45) ([**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) const & x) = default<br>_Construct the basis by copy._  |
|   | [**Impl**](#function-impl-55) ([**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) && x) = default<br>_Construct the basis from an r-value._  |
|  void | [**eval\_basis**](#function-eval_basis) (host\_t&lt; DField&lt; IdxRange&lt; DDim &gt; &gt; &gt; values, Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & x) const<br>_Evaluate the basis at the given coordinate._  |
|  [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) & | [**operator=**](#function-operator) ([**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) const & x) = default<br>_Copy-assign the class._  |
|  [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) & | [**operator=**](#function-operator_1) ([**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) && x) = default<br>_Move-assign the class._  |
|   | [**~Impl**](#function-impl) () = default<br> |




























## Detailed Description


The [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) class holds the implementation of the [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md).




**Template parameters:**


* `MemorySpace` Indicates where the object is saved. This is either on the host or the device. 




    
## Public Types Documentation




### typedef discrete\_dimension\_type 

_The tag which identifies the basis._ 
```C++
using TriangularBernsteinPolynomialBasis< X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D >::Impl< DDim, MemorySpace >::discrete_dimension_type =  TriangularBernsteinPolynomialBasis;
```




<hr>



### typedef discrete\_domain\_type 

_The type of the index range of the basis._ 
```C++
using TriangularBernsteinPolynomialBasis< X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D >::Impl< DDim, MemorySpace >::discrete_domain_type =  IdxRange<DDim>;
```




<hr>



### typedef discrete\_element\_type 

_The type of an index of an element of the basis._ 
```C++
using TriangularBernsteinPolynomialBasis< X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D >::Impl< DDim, MemorySpace >::discrete_element_type =  Idx<DDim>;
```




<hr>



### typedef discrete\_vector\_type 

_The type of an index step from one element of the basis to another._ 
```C++
using TriangularBernsteinPolynomialBasis< X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D >::Impl< DDim, MemorySpace >::discrete_vector_type =  IdxStep<DDim>;
```




<hr>
## Public Functions Documentation




### function Impl [2/5]

_Construct the basis from the barycentric coordinate mapping._ 
```C++
inline explicit TriangularBernsteinPolynomialBasis::Impl::Impl (
    CartesianToBarycentric < X , Y , Corner1Tag, Corner2Tag, Corner3Tag > const & coord_changer
) 
```





**Parameters:**


* `coord_changer` The class which converts Cartesian coordinates to barycentric coordinates. 




        

<hr>



### function Impl [3/5]

_Construct the basis by copy. This constructor is used to create the class on a different memory space._ 
```C++
template<class OriginMemorySpace>
inline explicit TriangularBernsteinPolynomialBasis::Impl::Impl (
    Impl < DDim, OriginMemorySpace > const & impl
) 
```





**Parameters:**


* `impl` The implementation of the origin memory space. 




        

<hr>



### function Impl [4/5]

_Construct the basis by copy._ 
```C++
TriangularBernsteinPolynomialBasis::Impl::Impl (
    Impl const & x
) = default
```





**Parameters:**


* `x` The basis to be copied. 




        

<hr>



### function Impl [5/5]

_Construct the basis from an r-value._ 
```C++
TriangularBernsteinPolynomialBasis::Impl::Impl (
    Impl && x
) = default
```





**Parameters:**


* `x` The temporary basis to be copied. 




        

<hr>



### function eval\_basis 

_Evaluate the basis at the given coordinate._ 
```C++
void TriangularBernsteinPolynomialBasis::Impl::eval_basis (
    host_t< DField< IdxRange< DDim > > > values,
    Coord< X , Y > const & x
) const
```





**Parameters:**


* `values` The values of the basis functions at the coordinate. 
* `x` The coordinate where the polynomials should be evaluated. 




        

<hr>



### function operator= 

_Copy-assign the class._ 
```C++
Impl & TriangularBernsteinPolynomialBasis::Impl::operator= (
    Impl const & x
) = default
```





**Parameters:**


* `x` An reference to another [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md). 



**Returns:**

A reference to this object. 





        

<hr>



### function operator= 

_Move-assign the class._ 
```C++
Impl & TriangularBernsteinPolynomialBasis::Impl::operator= (
    Impl && x
) = default
```





**Parameters:**


* `x` An rvalue to another [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md). 



**Returns:**

A reference to this object. 





        

<hr>



### function ~Impl 

```C++
TriangularBernsteinPolynomialBasis::Impl::~Impl () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/bernstein.hpp`

