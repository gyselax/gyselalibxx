

# Class InfNormResidual

**template &lt;class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**InfNormResidual**](classInfNormResidual.md)



_A Ginkgo stopping criterion based on the infinite norm of the residual._ [More...](#detailed-description)

* `#include <inf_norm_residual_criterion.hpp>`



Inherits the following classes: gko::EnablePolymorphicObject< InfNormResidual< Kokkos::DefaultExecutionSpace >, gko::stop::Criterion >


































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GKO\_CREATE\_FACTORY\_PARAMETERS**](#function-gko_create_factory_parameters) (parameters, Factory) <br>_The parameters which can be passed to the factory._  |
|   | [**GKO\_ENABLE\_BUILD\_METHOD**](#function-gko_enable_build_method) (Factory) <br>_Define the static build() method which returns the parameters used to create a Factory._  |
|   | [**GKO\_ENABLE\_CRITERION\_FACTORY**](#function-gko_enable_criterion_factory) ([**InfNormResidual**](classInfNormResidual.md)&lt; ExecSpace &gt;, parameters, Factory) <br>_Define the Factory class which creates the criterion, and the get\_parameters() accessor._  |
|  bool | [**check\_impl**](#function-check_impl) (gko::uint8 stopping\_id, bool set\_finalised, gko::array&lt; gko::stopping\_status &gt; \* stop\_status, bool \* one\_changed, gko::stop::Criterion::Updater const & updater) override<br>_Check if the infinite norm of the residual is below the threshold._  |
























## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**InfNormResidual**](#function-infnormresidual-12) (std::shared\_ptr&lt; gko::Executor const &gt; exec) <br>_Constructor used by gko::EnablePolymorphicObject._  |
|   | [**InfNormResidual**](#function-infnormresidual-22) (Factory const \* factory, gko::stop::CriterionArgs const & args) <br>_Constructor used by the factory. Computes the convergence threshold from the right-hand side._  |




## Detailed Description


The criterion is satisfied when \|\|r\|\|\_inf &lt;= reduction\_factor \* \|\|b\|\|\_inf, where r is the residual vector provided by the solver at each iteration (for CG/BiCGSTAB this is the residual computed by the recurrence of the method, not an explicitly re-computed b - Ax). Only single right-hand sides are supported.




**Template parameters:**


* `ExecSpace` The Kokkos execution space which has access to the Ginkgo executor memory. 




    
## Public Functions Documentation




### function GKO\_CREATE\_FACTORY\_PARAMETERS 

_The parameters which can be passed to the factory._ 
```C++
inline InfNormResidual::GKO_CREATE_FACTORY_PARAMETERS (
    parameters,
    Factory
) 
```




<hr>



### function GKO\_ENABLE\_BUILD\_METHOD 

_Define the static build() method which returns the parameters used to create a Factory._ 
```C++
InfNormResidual::GKO_ENABLE_BUILD_METHOD (
    Factory
) 
```




<hr>



### function GKO\_ENABLE\_CRITERION\_FACTORY 

_Define the Factory class which creates the criterion, and the get\_parameters() accessor._ 
```C++
InfNormResidual::GKO_ENABLE_CRITERION_FACTORY (
    InfNormResidual < ExecSpace >,
    parameters,
    Factory
) 
```




<hr>



### function check\_impl 

_Check if the infinite norm of the residual is below the threshold._ 
```C++
inline bool InfNormResidual::check_impl (
    gko::uint8 stopping_id,
    bool set_finalised,
    gko::array< gko::stopping_status > * stop_status,
    bool * one_changed,
    gko::stop::Criterion::Updater const & updater
) override
```



Note: This function should be private but it is public due to CUDA restrictions.




**Parameters:**


* `stopping_id` The id of the stopping criterion, saved in the stopping status. 
* `set_finalised` Controls if the current version should count as finalised. 
* `stop_status` The status of each right-hand side. 
* `one_changed` Indicates if the status of a right-hand side was modified. 
* `updater` The object containing the solver's current state.



**Returns:**

True if the criterion is satisfied. 





        

<hr>
## Protected Functions Documentation




### function InfNormResidual [1/2]

_Constructor used by gko::EnablePolymorphicObject._ 
```C++
inline explicit InfNormResidual::InfNormResidual (
    std::shared_ptr< gko::Executor const > exec
) 
```





**Parameters:**


* `exec` The Ginkgo executor. 




        

<hr>



### function InfNormResidual [2/2]

_Constructor used by the factory. Computes the convergence threshold from the right-hand side._ 
```C++
inline explicit InfNormResidual::InfNormResidual (
    Factory const * factory,
    gko::stop::CriterionArgs const & args
) 
```





**Parameters:**


* `factory` The factory holding the parameters. 
* `args` The arguments of the solver (system matrix, right-hand side, ...). 




        

<hr>## Friends Documentation





### friend EnablePolymorphicObject&lt; InfNormResidual&lt; ExecSpace &gt;, gko::stop::Criterion &gt; 

```C++
class InfNormResidual::EnablePolymorphicObject< InfNormResidual< ExecSpace >, gko::stop::Criterion > (
    gko::EnablePolymorphicObject< InfNormResidual < ExecSpace >, gko::stop::Criterion >
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/inf_norm_residual_criterion.hpp`

