

# Namespace ndEval



[**Namespace List**](namespaces.md) **>** [**ndEval**](namespacendEval.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**deriv**](#function-deriv) (Evaluator const & evaluator, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeEval, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutOut &gt; output, ConstField&lt; ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff &gt; coeffs) <br>_Compute first-order partial derivatives of a scalar function as a_ [_**VectorField**_](classVectorField.md) _._ |
|  void | [**deriv**](#function-deriv) (Evaluator const & evaluator, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeEval, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutOut &gt; output, ConstField&lt; CoordType, IdxRangeCoord, MemorySpace, LayoutCoords &gt; coords, ConstField&lt; ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff &gt; coeffs) <br>_Compute first-order partial derivatives at explicit coordinates as a_ [_**VectorField**_](classVectorField.md) _._ |
|  KOKKOS\_FUNCTION [**Vector**](classTensor.md)&lt; ElementType, VectorDims... &gt; | [**evaluate**](#function-evaluate) (Evaluator const & evaluator, CoordType const & coord, [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeCoeff, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutCoeff &gt; coeffs) <br>_Evaluate an ND interpolation at a single coordinate on each component of a_ [_**VectorField**_](classVectorField.md) _._ |
|  void | [**evaluate**](#function-evaluate) (Evaluator const & evaluator, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeEval, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutOut &gt; output, [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeCoeff, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutCoeff &gt; coeffs) <br>_Evaluate an ND interpolation on each component of a_ [_**VectorField**_](classVectorField.md) _._ |
|  void | [**evaluate**](#function-evaluate) (Evaluator const & evaluator, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeEval, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutOut &gt; output, ConstField&lt; CoordType, IdxRangeCoord, MemorySpace, LayoutCoords &gt; coords, [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeCoeff, VectorIndexSet&lt; VectorDims... &gt;, MemorySpace, LayoutCoeff &gt; coeffs) <br>_Evaluate an ND interpolation at explicit coordinates on each component of a_ [_**VectorField**_](classVectorField.md) _._ |




























## Public Functions Documentation




### function deriv 

_Compute first-order partial derivatives of a scalar function as a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<concepts::InterpolationEvaluator Evaluator, class ElementType, class IdxRangeEval, class IdxRangeCoeff, class MemorySpace, class LayoutOut, class LayoutCoeff, class... VectorDims>
void ndEval::deriv (
    Evaluator const & evaluator,
    VectorField < ElementType, IdxRangeEval, VectorIndexSet< VectorDims... >, MemorySpace, LayoutOut > output,
    ConstField< ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff > coeffs
) 
```



For each dimension `VectorDim` in `VectorDims`..., evaluates the first derivative of the scalar function defined by `coeffs` in the direction `VectorDim::Dual`, storing the result in the corresponding component of `output:` 
```C++
evaluator.deriv(Idx<ddc::Deriv<typename VectorDim::Dual>>(1),
                ddcHelper::get<VectorDim>(output), coeffs);
```



This is the [**VectorField**](classVectorField.md) counterpart of `evaluator.deriv(deriv_order, eval, coeffs)`, computing partial derivatives in all directions simultaneously.




**Template parameters:**


* `Evaluator` A class satisfying `concepts::InterpolationEvaluator`. 
* `VectorDims` The dimensions labelling the vector components. Each `VectorDim` must expose a `Dual` type identifying the continuous dimension in which the derivative is taken.



**Parameters:**


* `evaluator` The ND interpolation evaluator. 
* `output` On output, each component holds the first partial derivative of the scalar function in the direction `VectorDim::Dual`. 
* `coeffs` The scalar interpolation coefficients. 




        

<hr>



### function deriv 

_Compute first-order partial derivatives at explicit coordinates as a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<concepts::InterpolationEvaluator Evaluator, class ElementType, class IdxRangeEval, class CoordType, class IdxRangeCoord, class IdxRangeCoeff, class MemorySpace, class LayoutOut, class LayoutCoords, class LayoutCoeff, class... VectorDims>
void ndEval::deriv (
    Evaluator const & evaluator,
    VectorField < ElementType, IdxRangeEval, VectorIndexSet< VectorDims... >, MemorySpace, LayoutOut > output,
    ConstField< CoordType, IdxRangeCoord, MemorySpace, LayoutCoords > coords,
    ConstField< ElementType, IdxRangeCoeff, MemorySpace, LayoutCoeff > coeffs
) 
```



For each dimension `VectorDim` in `VectorDims`..., evaluates the first derivative of the scalar function defined by `coeffs` in the direction `VectorDim::Dual` at positions given by `coords`, storing the result in the corresponding component of `output:` 
```C++
evaluator.deriv(Idx<ddc::Deriv<typename VectorDim::Dual>>(1),
                ddcHelper::get<VectorDim>(output), coords, coeffs);
```



This is the [**VectorField**](classVectorField.md) counterpart of `evaluator.deriv(deriv_order, eval, coords, coeffs)`.




**Template parameters:**


* `Evaluator` A class satisfying `concepts::InterpolationEvaluator`. 
* `VectorDims` The dimensions labelling the vector components. Each `VectorDim` must expose a `Dual` type identifying the continuous dimension in which the derivative is taken.



**Parameters:**


* `evaluator` The ND interpolation evaluator. 
* `output` On output, each component holds the first partial derivative of the scalar function in the direction `VectorDim::Dual`. 
* `coords` The evaluation coordinates. 
* `coeffs` The scalar interpolation coefficients. 




        

<hr>



### function evaluate 

_Evaluate an ND interpolation at a single coordinate on each component of a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<concepts::InterpolationEvaluator Evaluator, class CoordType, class ElementType, class IdxRangeCoeff, class MemorySpace, class LayoutCoeff, class... VectorDims>
KOKKOS_FUNCTION Vector < ElementType, VectorDims... > ndEval::evaluate (
    Evaluator const & evaluator,
    CoordType const & coord,
    VectorField < const ElementType, IdxRangeCoeff, VectorIndexSet< VectorDims... >, MemorySpace, LayoutCoeff > coeffs
) 
```



For each dimension `VectorDim` in `VectorDims`..., calls 
```C++
evaluator(coord, ddcHelper::get<VectorDim>(coeffs));
```
 and packs the results into the returned Vector.


This is the [**VectorField**](classVectorField.md) counterpart of `evaluator.operator()(coord, coeffs)`.




**Template parameters:**


* `Evaluator` A class satisfying `concepts::InterpolationEvaluator`. 
* `VectorDims` The dimensions labelling the vector components (passed to `ddcHelper::get` to extract the corresponding scalar field).



**Parameters:**


* `evaluator` The ND interpolation evaluator. 
* `coord` The coordinate at which the interpolation is evaluated. 
* `coeffs` The interpolation coefficients for each vector component.



**Returns:**

The interpolated value of each vector component at `coord`. 





        

<hr>



### function evaluate 

_Evaluate an ND interpolation on each component of a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<concepts::InterpolationEvaluator Evaluator, class ElementType, class IdxRangeEval, class IdxRangeCoeff, class MemorySpace, class LayoutOut, class LayoutCoeff, class... VectorDims>
void ndEval::evaluate (
    Evaluator const & evaluator,
    VectorField < ElementType, IdxRangeEval, VectorIndexSet< VectorDims... >, MemorySpace, LayoutOut > output,
    VectorField < const ElementType, IdxRangeCoeff, VectorIndexSet< VectorDims... >, MemorySpace, LayoutCoeff > coeffs
) 
```



For each dimension `VectorDim` in `VectorDims`..., calls 
```C++
evaluator(ddcHelper::get<VectorDim>(output), ddcHelper::get<VectorDim>(coeffs));
```



This is the [**VectorField**](classVectorField.md) counterpart of `evaluator.operator()(eval, coeffs)`.




**Template parameters:**


* `Evaluator` A class satisfying `concepts::InterpolationEvaluator`. 
* `VectorDims` The dimensions labelling the vector components (passed to `ddcHelper::get` to extract the corresponding scalar field).



**Parameters:**


* `evaluator` The ND interpolation evaluator. 
* `output` On output, the interpolated values for each vector component. 
* `coeffs` The interpolation coefficients for each vector component. 




        

<hr>



### function evaluate 

_Evaluate an ND interpolation at explicit coordinates on each component of a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<concepts::InterpolationEvaluator Evaluator, class ElementType, class IdxRangeEval, class CoordType, class IdxRangeCoord, class IdxRangeCoeff, class MemorySpace, class LayoutOut, class LayoutCoords, class LayoutCoeff, class... VectorDims>
void ndEval::evaluate (
    Evaluator const & evaluator,
    VectorField < ElementType, IdxRangeEval, VectorIndexSet< VectorDims... >, MemorySpace, LayoutOut > output,
    ConstField< CoordType, IdxRangeCoord, MemorySpace, LayoutCoords > coords,
    VectorField < const ElementType, IdxRangeCoeff, VectorIndexSet< VectorDims... >, MemorySpace, LayoutCoeff > coeffs
) 
```



For each dimension `VectorDim` in `VectorDims`..., calls 
```C++
evaluator(ddcHelper::get<VectorDim>(output), coords, ddcHelper::get<VectorDim>(coeffs));
```



This is the [**VectorField**](classVectorField.md) counterpart of `evaluator.operator()(eval, coords, coeffs)`.




**Template parameters:**


* `Evaluator` A class satisfying `concepts::InterpolationEvaluator`. 
* `VectorDims` The dimensions labelling the vector components (passed to `ddcHelper::get` to extract the corresponding scalar field).



**Parameters:**


* `evaluator` The ND interpolation evaluator. 
* `output` On output, the interpolated values for each vector component. 
* `coords` The evaluation coordinates. 
* `coeffs` The interpolation coefficients for each vector component. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/vector_field_evaluation.hpp`

