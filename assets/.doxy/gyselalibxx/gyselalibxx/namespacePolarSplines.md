

# Namespace PolarSplines



[**Namespace List**](namespaces.md) **>** [**PolarSplines**](namespacePolarSplines.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  DField&lt; IdxRange&lt; DDim &gt;, MemorySpace &gt; | [**integrals**](#function-integrals) (ExecSpace const & execution\_space, DField&lt; IdxRange&lt; DDim &gt;, MemorySpace &gt; int\_vals) <br>_Get the integrals over the logical domain of each of the polar B-splines._  |




























## Public Functions Documentation




### function integrals 

_Get the integrals over the logical domain of each of the polar B-splines._ 
```C++
template<class ExecSpace, class DDim, class MemorySpace>
DField< IdxRange< DDim >, MemorySpace > PolarSplines::integrals (
    ExecSpace const & execution_space,
    DField< IdxRange< DDim >, MemorySpace > int_vals
) 
```





**Parameters:**


* `execution_space` The execution space on which the integrals should be calculated. 
* `int_vals` The values of the integrals over the logical domain of each of the polar B-splines. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_bsplines.hpp`

