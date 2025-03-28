

# Class CartesianToCzarny

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class [**R**](structR.md), class [**Theta**](structTheta.md)&gt;**



[**ClassList**](annotated.md) **>** [**CartesianToCzarny**](classCartesianToCzarny.md)



_A class for describing the Czarny 2D mapping._ [More...](#detailed-description)

* `#include <cartesian_to_czarny.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef [**X**](structX.md) | [**cartesian\_tag\_x**](#typedef-cartesian_tag_x)  <br>_Indicate the first physical coordinate._  |
| typedef [**Y**](structY.md) | [**cartesian\_tag\_y**](#typedef-cartesian_tag_y)  <br>_Indicate the second physical coordinate._  |
| typedef [**R**](structR.md) | [**curvilinear\_tag\_r**](#typedef-curvilinear_tag_r)  <br>_Indicate the first logical coordinate._  |
| typedef [**Theta**](structTheta.md) | [**curvilinear\_tag\_theta**](#typedef-curvilinear_tag_theta)  <br>_Indicate the second logical coordinate._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CartesianToCzarny**](#function-cartesiantoczarny-13) (double epsilon, double e) <br>_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from parameters._ |
|  KOKKOS\_FUNCTION | [**CartesianToCzarny**](#function-cartesiantoczarny-23) ([**CartesianToCzarny**](classCartesianToCzarny.md) const & other) <br>_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(lvalue)._ |
|   | [**CartesianToCzarny**](#function-cartesiantoczarny-33) ([**CartesianToCzarny**](classCartesianToCzarny.md) && x) = default<br>_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another temporary_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(rvalue)._ |
|  KOKKOS\_FUNCTION double | [**e**](#function-e) () const<br>_Return the_  _parameter._ |
|  KOKKOS\_FUNCTION double | [**epsilon**](#function-epsilon) () const<br>_Return the_  _parameter._ |
|  [**CzarnyToCartesian**](classCzarnyToCartesian.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md), [**X**](structX.md), [**Y**](structY.md) &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**operator()**](#function-operator) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & coord) const<br>_Convert the coordinate (x,y) to the equivalent_  _coordinate._ |
|  [**CartesianToCzarny**](classCartesianToCzarny.md) & | [**operator=**](#function-operator_1) ([**CartesianToCzarny**](classCartesianToCzarny.md) const & x) = default<br>_Assign a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(lvalue)._ |
|  [**CartesianToCzarny**](classCartesianToCzarny.md) & | [**operator=**](#function-operator_2) ([**CartesianToCzarny**](classCartesianToCzarny.md) && x) = default<br>_Assign a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another temporary_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(rvalue)._ |
|   | [**~CartesianToCzarny**](#function-cartesiantoczarny) () = default<br> |




























## Detailed Description


The mapping  is defined by








with  and  and  given as parameters. 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::CoordArg =  Coord<X, Y>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::CoordResult =  Coord<R, Theta>;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first physical coordinate._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second physical coordinate._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::cartesian_tag_y =  Y;
```




<hr>



### typedef curvilinear\_tag\_r 

_Indicate the first logical coordinate._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::curvilinear_tag_r =  R;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using CartesianToCzarny< X, Y, R, Theta >::curvilinear_tag_theta =  Theta;
```




<hr>
## Public Functions Documentation




### function CartesianToCzarny [1/3]

_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from parameters._
```C++
inline CartesianToCzarny::CartesianToCzarny (
    double epsilon,
    double e
) 
```





**Parameters:**


* `epsilon` The  parameter in the definition of the mapping [**CartesianToCzarny**](classCartesianToCzarny.md).
* `e` The  parameter in the definition of the mapping [**CartesianToCzarny**](classCartesianToCzarny.md).



**See also:** [**CartesianToCzarny**](classCartesianToCzarny.md) 



        

<hr>



### function CartesianToCzarny [2/3]

_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(lvalue)._
```C++
inline KOKKOS_FUNCTION CartesianToCzarny::CartesianToCzarny (
    CartesianToCzarny const & other
) 
```





**Parameters:**


* `other` [**CartesianToCzarny**](classCartesianToCzarny.md) mapping used to instantiate the new one. 




        

<hr>



### function CartesianToCzarny [3/3]

_Instantiate a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another temporary_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(rvalue)._
```C++
CartesianToCzarny::CartesianToCzarny (
    CartesianToCzarny && x
) = default
```





**Parameters:**


* `x` [**CartesianToCzarny**](classCartesianToCzarny.md) mapping used to instantiate the new one. 




        

<hr>



### function e 

_Return the_  _parameter._
```C++
inline KOKKOS_FUNCTION double CartesianToCzarny::e () const
```





**Returns:**

The value of .




**See also:** [**CartesianToCzarny**](classCartesianToCzarny.md) 



        

<hr>



### function epsilon 

_Return the_  _parameter._
```C++
inline KOKKOS_FUNCTION double CartesianToCzarny::epsilon () const
```





**Returns:**

The value of .




**See also:** [**CartesianToCzarny**](classCartesianToCzarny.md) 



        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CzarnyToCartesian < R , Theta , X , Y > CartesianToCzarny::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function operator() 

_Convert the coordinate (x,y) to the equivalent_  _coordinate._
```C++
inline KOKKOS_FUNCTION Coord< R , Theta > CartesianToCzarny::operator() (
    Coord< X , Y > const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(lvalue)._
```C++
CartesianToCzarny & CartesianToCzarny::operator= (
    CartesianToCzarny const & x
) = default
```





**Parameters:**


* `x` [**CartesianToCzarny**](classCartesianToCzarny.md) mapping used to assign.



**Returns:**

The [**CartesianToCzarny**](classCartesianToCzarny.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCzarny**_](classCartesianToCzarny.md) _from another temporary_[_**CartesianToCzarny**_](classCartesianToCzarny.md) _(rvalue)._
```C++
CartesianToCzarny & CartesianToCzarny::operator= (
    CartesianToCzarny && x
) = default
```





**Parameters:**


* `x` [**CartesianToCzarny**](classCartesianToCzarny.md) mapping used to assign.



**Returns:**

The [**CartesianToCzarny**](classCartesianToCzarny.md) assigned. 





        

<hr>



### function ~CartesianToCzarny 

```C++
CartesianToCzarny::~CartesianToCzarny () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/cartesian_to_czarny.hpp`

