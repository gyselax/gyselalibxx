

# Class Matrix\_Corner\_Block



[**ClassList**](annotated.md) **>** [**Matrix\_Corner\_Block**](classMatrix__Corner__Block.md)



_A class representing a matrix with the following block pattern:_ [More...](#detailed-description)

* `#include <matrix_corner_block.hpp>`



Inherits the following classes: [Matrix](classMatrix.md)


Inherited by the following classes: [Matrix\_Centre\_Block](classMatrix__Centre__Block.md),  [Matrix\_Periodic\_Banded](classMatrix__Periodic__Banded.md)




















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Corner\_Block**](#function-matrix_corner_block-12) (int n, int k, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q) <br>_A constructor for the matrix._  |
| virtual void | [**factorise**](#function-factorise) () override<br>_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._  |
| virtual double | [**get\_element**](#function-get_element) (int i, int j) override const<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](#function-set_element) (int i, int j, double a\_ij) override<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](#function-solve_inplace) (DSpan1D b) override const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](#function-solve_multiple_inplace) (DSpan2D bx) override const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](#function-solve_transpose_inplace) (DSpan1D b) override const<br>_Solve the transposed matrix equation in place._  |


## Public Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|   | [**Matrix**](classMatrix.md#function-matrix) (int mat\_size) <br>_A constructor for the matrix._  |
| virtual void | [**factorise**](classMatrix.md#function-factorise) () <br>_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._  |
| virtual double | [**get\_element**](classMatrix.md#function-get_element) (int i, int j) const = 0<br>_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._ |
|  int | [**get\_size**](classMatrix.md#function-get_size) () const<br>_Get the size of the n x n_ [_**Matrix**_](classMatrix.md) _._ |
| virtual void | [**set\_element**](classMatrix.md#function-set_element) (int i, int j, double a\_ij) = 0<br>_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._ |
| virtual DSpan1D | [**solve\_inplace**](classMatrix.md#function-solve_inplace) (DSpan1D b) const<br>_Solve the matrix equation in place._  |
| virtual DSpan2D | [**solve\_multiple\_inplace**](classMatrix.md#function-solve_multiple_inplace) (DSpan2D bx) const<br>_Solve multiple matrix equations in place._  |
| virtual DSpan1D | [**solve\_transpose\_inplace**](classMatrix.md#function-solve_transpose_inplace) (DSpan1D b) const<br>_Solve the transposed matrix equation in place._  |
| virtual  | [**~Matrix**](classMatrix.md#function-matrix) () = default<br> |




## Public Static Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_banded**](classMatrix.md#function-make_new_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a banded matrix. This method returns the appropriate subclass to minimise memory and computations._ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_block\_with\_banded\_region**](classMatrix.md#function-make_new_block_with_banded_region) (int n, int kl, int ku, bool pds, int block1\_size, int block2\_size=0) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This method returns the appropriate subclass to minimise memory and computations. This matrix must be able to be described by the following block matrices:_ |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**make\_new\_periodic\_banded**](classMatrix.md#function-make_new_periodic_banded) (int n, int kl, int ku, bool pds) <br>_Get a new_ [_**Matrix**_](classMatrix.md) _representing a periodic banded matrix. This method returns the appropriate subclass to minimise memory and computations. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_ |










## Protected Attributes

| Type | Name |
| ---: | :--- |
|  DSpan2D | [**Abm\_1\_gamma**](#variable-abm_1_gamma)  <br>_The upper-right sub-matrix. The sub-matrix gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. The sub-matrix beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._ |
|  std::unique\_ptr&lt; double[]&gt; | [**Abm\_1\_gamma\_ptr**](#variable-abm_1_gamma_ptr)  <br>_Data storage for the upper-right sub-matrix. This memory describes gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. This memory describes beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._ |
|  [**Matrix\_Dense**](classMatrix__Dense.md) | [**delta**](#variable-delta)  <br>_The sub-matrix delta which is a dense matrix._  |
|  int const | [**k**](#variable-k)  <br>_The size of the dense corner matrix delta._  |
|  DSpan2D | [**lambda**](#variable-lambda)  <br>_The sub-matrix lambda._  |
|  std::unique\_ptr&lt; double[]&gt; | [**lambda\_ptr**](#variable-lambda_ptr)  <br>_Data storage for the sub-matrix lambda._  |
|  int const | [**nb**](#variable-nb)  <br>_The size of the banded matrix._  |
|  std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; | [**q\_block**](#variable-q_block)  <br>_The sub-matrix Q which is a banded matrix._  |


## Protected Attributes inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
|  int const | [**n**](classMatrix.md#variable-n)  <br>_The matrix size._  |






























## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**Matrix\_Corner\_Block**](#function-matrix_corner_block-22) (int n, int k, std::unique\_ptr&lt; [**Matrix**](classMatrix.md) &gt; q, int lambda\_size1, int lambda\_size2) <br>_A constructor for the matrix._  |
| virtual void | [**calculate\_delta\_to\_factorise**](#function-calculate_delta_to_factorise) () <br>_Calculate the contents of the dense matrix_  _that will be factorised. This is the element_ _of the blockwise LU decomposition of the matrix._ |
| virtual DSpan1D | [**solve\_gamma\_section**](#function-solve_gamma_section) (DSpan1D const u, DView1D const v) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_gamma\_section\_transpose**](#function-solve_gamma_section_transpose) (DSpan1D const v, DView1D const u) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_lambda\_section**](#function-solve_lambda_section) (DSpan1D v, DView1D u) const<br>_Calculate the solution to the following equation:_  |
| virtual DSpan1D | [**solve\_lambda\_section\_transpose**](#function-solve_lambda_section_transpose) (DSpan1D u, DView1D v) const<br>_Calculate the solution to the following equation:_  |


## Protected Functions inherited from Matrix

See [Matrix](classMatrix.md)

| Type | Name |
| ---: | :--- |
| virtual int | [**factorise\_method**](classMatrix.md#function-factorise_method) () = 0<br>_Call the factorisation method._  |
| virtual int | [**solve\_inplace\_method**](classMatrix.md#function-solve_inplace_method) (double \* b, char transpose, int n\_equations) const = 0<br>_Call the LAPACK solve method._  |






## Detailed Description


 where Q is a banded matrix, and Q and delta are square matrices.


The matrix equation is factorised using a LU decomposition. The equation is then solved using this decomposition as described in section 2.5.2.1 of Emily Bourne's thesis [1].


[1] Emily Bourne, "Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code". December 2022. 


    
## Public Functions Documentation




### function Matrix\_Corner\_Block [1/2]

_A constructor for the matrix._ 
```C++
Matrix_Corner_Block::Matrix_Corner_Block (
    int n,
    int k,
    std::unique_ptr< Matrix > q
) 
```





**Parameters:**


* `n` The size of the n x n matrix. 
* `k` The size of the k x k sub-matrix delta. 
* `q` The sub-matrix Q. 




        

<hr>



### function factorise 

_Factorise the matrix. This method prepares the matrix for a call to the solve methods. For most matrix types a call to factorise causes the LU decomposition to be calculated. The elements of the matrix should not be accessed once this method has been called._ 
```C++
virtual void Matrix_Corner_Block::factorise () override
```



Implements [*Matrix::factorise*](classMatrix.md#function-factorise)


<hr>



### function get\_element 

_A method to get an element from the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual double Matrix_Corner_Block::get_element (
    int i,
    int j
) override const
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 



**Returns:**

The element at position (i,j). 





        
Implements [*Matrix::get\_element*](classMatrix.md#function-get_element)


<hr>



### function set\_element 

_A method to set an element of the_ [_**Matrix**_](classMatrix.md) _._
```C++
virtual void Matrix_Corner_Block::set_element (
    int i,
    int j,
    double a_ij
) override
```





**Parameters:**


* `i` The row index. 
* `j` The column index. 
* `a_ij` The value that the element should be set to. 




        
Implements [*Matrix::set\_element*](classMatrix.md#function-set_element)


<hr>



### function solve\_inplace 

_Solve the matrix equation in place._ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_inplace (
    DSpan1D b
) override const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix::solve\_inplace*](classMatrix.md#function-solve_inplace)


<hr>



### function solve\_multiple\_inplace 

_Solve multiple matrix equations in place._ 
```C++
virtual DSpan2D Matrix_Corner_Block::solve_multiple_inplace (
    DSpan2D bx
) override const
```



Solve the following matrix equation:  for multiple values of  and . The first dimension is iterated over with each slice representing an equation to be solved. The result  is saved into the memory allocated for .




**Parameters:**


* `bx` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix::solve\_multiple\_inplace*](classMatrix.md#function-solve_multiple_inplace)


<hr>



### function solve\_transpose\_inplace 

_Solve the transposed matrix equation in place._ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_transpose_inplace (
    DSpan1D b
) override const
```



Solve the following matrix equation:  The result  is saved into the memory allocated for .




**Parameters:**


* `b` The right-hand side of the equation on input. The solution on output. 



**Returns:**

The solution . 





        
Implements [*Matrix::solve\_transpose\_inplace*](classMatrix.md#function-solve_transpose_inplace)


<hr>
## Protected Attributes Documentation




### variable Abm\_1\_gamma 

_The upper-right sub-matrix. The sub-matrix gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. The sub-matrix beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._
```C++
DSpan2D Matrix_Corner_Block::Abm_1_gamma;
```




<hr>



### variable Abm\_1\_gamma\_ptr 

_Data storage for the upper-right sub-matrix. This memory describes gamma before the_ [_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called. This memory describes beta after the_[_**factorise()**_](classMatrix__Corner__Block.md#function-factorise) _method is called._
```C++
std::unique_ptr<double[]> Matrix_Corner_Block::Abm_1_gamma_ptr;
```




<hr>



### variable delta 

_The sub-matrix delta which is a dense matrix._ 
```C++
Matrix_Dense Matrix_Corner_Block::delta;
```




<hr>



### variable k 

_The size of the dense corner matrix delta._ 
```C++
int const Matrix_Corner_Block::k;
```




<hr>



### variable lambda 

_The sub-matrix lambda._ 
```C++
DSpan2D Matrix_Corner_Block::lambda;
```




<hr>



### variable lambda\_ptr 

_Data storage for the sub-matrix lambda._ 
```C++
std::unique_ptr<double[]> Matrix_Corner_Block::lambda_ptr;
```




<hr>



### variable nb 

_The size of the banded matrix._ 
```C++
int const Matrix_Corner_Block::nb;
```




<hr>



### variable q\_block 

_The sub-matrix Q which is a banded matrix._ 
```C++
std::unique_ptr<Matrix> Matrix_Corner_Block::q_block;
```




<hr>
## Protected Functions Documentation




### function Matrix\_Corner\_Block [2/2]

_A constructor for the matrix._ 
```C++
Matrix_Corner_Block::Matrix_Corner_Block (
    int n,
    int k,
    std::unique_ptr< Matrix > q,
    int lambda_size1,
    int lambda_size2
) 
```





**Parameters:**


* `n` The size of the n x n matrix. 
* `k` The size of the k x k sub-matrix delta. 
* `q` The sub-matrix Q. 
* `lambda_size1` The number of rows necessary to store the sub-matrix lambda. 
* `lambda_size2` The number of columns necessary to store the sub-matrix lambda. 




        

<hr>



### function calculate\_delta\_to\_factorise 

_Calculate the contents of the dense matrix_  _that will be factorised. This is the element_ _of the blockwise LU decomposition of the matrix._
```C++
virtual void Matrix_Corner_Block::calculate_delta_to_factorise () 
```



 It is defined as:





With  the solution of the equation:





which should already be stored in the variable Abm\_1\_gamma when this method is called. 


        

<hr>



### function solve\_gamma\_section 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_gamma_section (
    DSpan1D const u,
    DView1D const v
) const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the matrix equation. 
* `v` A vector of length k. This is an intermediate result of the calculation which solves the matrix equation.



**Returns:**

The vector u solution of the equation. 





        

<hr>



### function solve\_gamma\_section\_transpose 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_gamma_section_transpose (
    DSpan1D const v,
    DView1D const u
) const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `v` A vector of length k describing the final values in the right-hand side of the matrix to be solved. 
* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the transposed matrix equation.



**Returns:**

The vector v solution of the equation. 





        

<hr>



### function solve\_lambda\_section 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_lambda_section (
    DSpan1D v,
    DView1D u
) const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `v` A vector of length k describing the final values in the right-hand side of the matrix to be solved. 
* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the matrix equation.



**Returns:**

The vector v solution of the equation. 





        

<hr>



### function solve\_lambda\_section\_transpose 

_Calculate the solution to the following equation:_ 
```C++
virtual DSpan1D Matrix_Corner_Block::solve_lambda_section_transpose (
    DSpan1D u,
    DView1D v
) const
```






The result is saved into the argument . This calculation is one of the steps necessary to solve the matrix equation.




**Parameters:**


* `u` A vector of length (n-k). This is an intermediate result of the calculation which solves the transposed matrix equation. 
* `v` A vector of length k. This is an intermediate result of the calculation which solves the transposed matrix equation.



**Returns:**

The vector u solution of the equation. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_corner_block.hpp`

