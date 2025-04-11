

# Class AdvectionFieldFinder

**template &lt;class Mapping&gt;**



[**ClassList**](annotated.md) **>** [**AdvectionFieldFinder**](classAdvectionFieldFinder.md)



_Solve the Poisson-like equation and return the electric field for the coupled Vlasov equation._ [More...](#detailed-description)

* `#include <advection_field_rtheta.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::array&lt; std::array&lt; double, 2 &gt;, 2 &gt; | [**Matrix\_2x2**](#typedef-matrix_2x2)  <br>_Define a 2x2 matrix with an 2D array of an 2D array._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**AdvectionFieldFinder**](#function-advectionfieldfinder) (Mapping const & mapping, double const epsilon=1e-12) <br>_Instantiate a AdvectionFieldRTheta ._  |
|  void | [**operator()**](#function-operator) (host\_t&lt; DFieldRTheta &gt; electrostatic\_potential, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; advection\_field\_xy) const<br>_Compute the advection field from a Field of_  _values._ |
|  void | [**operator()**](#function-operator_1) (host\_t&lt; Spline2D &gt; electrostatic\_potential\_coef, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; advection\_field\_xy) const<br>_Compute the advection field from a spline representation of_  _solution. The B-splines basis used is the cross-product of two 1D B-splines basis._ |
|  void | [**operator()**](#function-operator_2) (host\_t&lt; [**PolarSplineMemRTheta**](structPolarSplineMem.md) &gt; & electrostatic\_potential\_coef, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; advection\_field\_xy) const<br>_Compute the advection field from the Poisson-like equation solution. The B-splines basis used is the polar B-splines (_ [_**PolarSplineMem**_](structPolarSplineMem.md) _)._ |
|  void | [**operator()**](#function-operator_3) (host\_t&lt; DFieldRTheta &gt; electrostatic\_potential, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; &gt; advection\_field\_rtheta, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; & advection\_field\_xy\_centre) const<br>_Compute the advection field from a Field of_  _values._ |
|  void | [**operator()**](#function-operator_4) (host\_t&lt; Spline2D &gt; electrostatic\_potential\_coef, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; &gt; advection\_field\_rtheta, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; & advection\_field\_xy\_centre) const<br>_Compute the advection field from a spline representation of_  _. The B-splines basis used is the cross-product of two 1D B-splines basis._ |
|  void | [**operator()**](#function-operator_5) (host\_t&lt; [**PolarSplineMemRTheta**](structPolarSplineMem.md) &gt; & electrostatic\_potential\_coef, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; &gt; advection\_field\_rtheta, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; & advection\_field\_xy\_centre) const<br>_Compute the advection field from the Poisson-like equation. The B-splines basis used is the polar B-splines (_ [_**PolarSplineMem**_](structPolarSplineMem.md) _)._ |




























## Detailed Description


The Vlasov-Poisson equations are given by



* (1) ,
* (2) ,
* (3) .




The functions are defined on a logical domain, and the mapping from the logical domain to the physical domain is written .


We here focus on equation (3). The  is already computed on B-splines with the given Poisson solver. Then in the AdvectionFieldFinder::operator() we compute the advection field ( ) thanks to (3) using the B-splines coefficients. Depending on the given mapping, the computation at the centre point is not always well-defined so we linearise around the centre point as explained in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).


The advection field can be computed along the logical domain axis or the physical domain axis.


1- In the first case, we compute the electric field thanks to (3) and
* ,
* ,




For ,  is ill-defined so we linearise  ,


with  computed thanks to



* ,
* ,




where  and  correspond to linearly independent directions.


Then the advection field along the physical domain axis is given by .


2- In the second case, the advection field along the logical domain axis is computed with
* ,
* with , the coefficients of the inverse metric tensor,
* , the coefficients of the metric tensor,
* , the unnormalized local covariants vectors.




Then, we compute  and , with  the Jacobian matrix of the transformation .


The equation (1) is solved thanks to advection operator (IAdvectionRTheta).




**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates.



**See also:** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) 



    
## Public Types Documentation




### typedef Matrix\_2x2 

_Define a 2x2 matrix with an 2D array of an 2D array._ 
```C++
using AdvectionFieldFinder< Mapping >::Matrix_2x2 =  std::array<std::array<double, 2>, 2>;
```




<hr>
## Public Functions Documentation




### function AdvectionFieldFinder 

_Instantiate a AdvectionFieldRTheta ._ 
```C++
inline explicit AdvectionFieldFinder::AdvectionFieldFinder (
    Mapping const & mapping,
    double const epsilon=1e-12
) 
```





**Parameters:**


* `mapping` The mapping  from the logical domain to the physical domain. 
* `epsilon` The parameter  for the linearisation of the electric field. 




        

<hr>



### function operator() 

_Compute the advection field from a Field of_  _values._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< DFieldRTheta > electrostatic_potential,
    host_t< DVectorFieldRTheta < X , Y > > advection_field_xy
) const
```





**Parameters:**


* `electrostatic_potential` The values of the solution  of the Poisson-like equation (2). 
* `advection_field_xy` The advection field on the physical axis. 




        

<hr>



### function operator() 

_Compute the advection field from a spline representation of_  _solution. The B-splines basis used is the cross-product of two 1D B-splines basis._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< Spline2D > electrostatic_potential_coef,
    host_t< DVectorFieldRTheta < X , Y > > advection_field_xy
) const
```





**Parameters:**


* `electrostatic_potential_coef` The spline representation of the solution  of the Poisson-like equation (2). 
* `advection_field_xy` The advection field on the physical axis. 




        

<hr>



### function operator() 

_Compute the advection field from the Poisson-like equation solution. The B-splines basis used is the polar B-splines (_ [_**PolarSplineMem**_](structPolarSplineMem.md) _)._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< PolarSplineMemRTheta > & electrostatic_potential_coef,
    host_t< DVectorFieldRTheta < X , Y > > advection_field_xy
) const
```





**Parameters:**


* `electrostatic_potential_coef` The polar spline representation of the solution  of the Poisson-like equation (2). 
* `advection_field_xy` The advection field on the physical axis. 




        

<hr>



### function operator() 

_Compute the advection field from a Field of_  _values._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< DFieldRTheta > electrostatic_potential,
    host_t< DVectorFieldRTheta < R , Theta > > advection_field_rtheta,
    DVector < X , Y > & advection_field_xy_centre
) const
```





**Parameters:**


* `electrostatic_potential` The values of the solution  of the Poisson-like equation (2). 
* `advection_field_rtheta` The advection field on the logical axis. It is expressed on the contravariant basis. 
* `advection_field_xy_centre` The advection field at the centre point on the Cartesian basis. 




        

<hr>



### function operator() 

_Compute the advection field from a spline representation of_  _. The B-splines basis used is the cross-product of two 1D B-splines basis._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< Spline2D > electrostatic_potential_coef,
    host_t< DVectorFieldRTheta < R , Theta > > advection_field_rtheta,
    DVector < X , Y > & advection_field_xy_centre
) const
```





**Parameters:**


* `electrostatic_potential_coef` The spline representation of the solution  of the Poisson-like equation (2). 
* `advection_field_rtheta` The advection field on the logical axis. It is expressed on the contravariant basis. 
* `advection_field_xy_centre` The advection field at the centre point on the Cartesian basis. 




        

<hr>



### function operator() 

_Compute the advection field from the Poisson-like equation. The B-splines basis used is the polar B-splines (_ [_**PolarSplineMem**_](structPolarSplineMem.md) _)._
```C++
inline void AdvectionFieldFinder::operator() (
    host_t< PolarSplineMemRTheta > & electrostatic_potential_coef,
    host_t< DVectorFieldRTheta < R , Theta > > advection_field_rtheta,
    DVector < X , Y > & advection_field_xy_centre
) const
```





**Parameters:**


* `electrostatic_potential_coef` The polar spline representation of the solution  of the Poisson-like equation (2). 
* `advection_field_rtheta` The advection field on the logical axis. It is expressed on the contravariant basis. 
* `advection_field_xy_centre` The advection field at the centre point on the Cartesian basis. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/advection_field/advection_field_rtheta.hpp`

