# Common compilation problems

<a name="KOKKOS-LAMBDA"></a>
## The closure type for a lambda cannot be used in the template argument type of a '\_\_global\_\_' function

When using `ddc::parallel_for_each` or `ddc::parallel_transform_reduce` the last argument passed to the functions is usually a lambda function. The usual syntax for lambda functions is:
```cpp
[capture_expression](ArgTypes arguments) { code; }
```
When compiling for GPU this syntax must be modified slightly.

The only valid capture expression is capture by copy (i.e. everything used in the code which is defined elsewhere is passed by copy). This is so that objects (especially scalars) can be copied to the device when they are needed. A capture by reference would lead to the code trying to use objects on the GPU which are only accessible from the CPU.

### Error Message

When this rule is not followed you will see the following error message:
```
${GYSELALIBXX_HOME}/vendor/kokkos/core/src/Cuda/Kokkos_Cuda_KernelLaunch.hpp(345): error: The closure type for a lambda ("lambda [](ArgTypes)->void", defined at <PATH_TO_FILE>:<LINE_NUMBER>) cannot be used in the template argument type of a __global__ function template instantiation, unless the lambda is defined within a __device__ or __global__ function, or the flag '-extended-lambda' is specified and the lambda is an extended lambda (a __device__ or __host__ __device__ lambda defined within a __host__ or __host__ __device__ function)
```

## Implicit capture of 'this' in extended lambda expression

Lambdas capture variables from their environment in order to use them in the expressions that are executed in the lambda function. They can capture local variables directly (by copy or reference), however in order to capture class variables the entire class must be captured. As captures for GPU require copies (see [above](#KOKKOS-LAMBDA)) such a copy could be extremely costly. It is therefore advised to avoid class capture.

Luckily it is simple to avoid the class being captured. It is sufficient to define local copies or accessors (i.e. `ddc::ChunkSpan`) to the variables that will be used in the function.

E.g. for the following code:
```cpp
double my_function() const
{
    return ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            m_class_field_x.domain(),
            KOKKOS_LAMBDA(IdxX ix) { return m_class_field_x(ix); });
}
```
where `m_class_field_x` is a class member function the code should be:
```cpp
double my_function() const
{
    DSpanX class_field_x_proxy = m_class_field_x.span_view();
    return ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            class_field_x_proxy.domain(),
            KOKKOS_LAMBDA(IdxX ix) { return class_field_x_proxy(ix); });
}
```

Unfortunately the error message is only a warning so if this is missed the code will still compile but the results may be wrong.

### Error Message

When this rule is not followed you will see the following error message:
```
<PATH_TO_FILE>(<LINE_NUMBER>): warning #20178-D: Implicit capture of 'this' in extended lambda expression
```

## Accessing allocated data

DDC provides 3 types of objects for storing/accessing array data:
-  `ddc::Chunk`
-  `ddc::ChunkSpan`
-  `ddc::ChunkView`

Whenever a new `ddc::Chunk` is created, memory is allocated. This includes when a `ddc::Chunk` is copied. In order to avoid accidental memory allocation the copy operator of `ddc::Chunk` is therefore deleted.

This is a little awkward as it means that a `ddc::Chunk` then cannot be used as a function argument or used inside a GPU function (as captures for GPU require copies, see [above](#KOKKOS-LAMBDA)). In order to get round this restriction the `ddc::ChunkSpan` type must be used. This type provides all the functionalities to access the data that is stored in a `ddc::Chunk` but creating/copying this object never allocates memory so there is no risk of accidental memory allocation. A `ddc::ChunkView` is almost identical to a `ddc::ChunkSpan` but the data it refers to cannot be modified.

In practice this means that while `ddc::Chunk` is used to initialise the memory, only `ddc::ChunkSpan` or `ddc::ChunkView` should be used to interact with the memory.

### Error Message

When this rule is not followed you may see a large number of errors. The following is a selection of the error messages that you may see:

The first error to appear is not very explicit. It simply indicates that the lvalue (the value on the left hand side of the assignment) has an unknown type. As a result the compiler does not know how to assign a value to it:
```
<PATH_TO_FILE>(<LINE_NUMBER>): error: expression must be a modifiable lvalue
```

The next error highlights the problem. It says that a `ddc::Chunk` cannot be copied using a constructor. This is because such an operation would allocate memory and this is not what is wanted.
```
<PATH_TO_FILE>(<LINE_NUMBER>): error: function "ddc::Chunk<ElementType, ddc::DiscreteDomain<DDims...>, Allocator>::Chunk(const ddc::Chunk<ElementType, ddc::DiscreteDomain<DDims...>, Allocator> &) [with ElementType=double, DDims=<Dim1, Dim2,...>, Allocator=ddc::KokkosAllocator<double, Kokkos::CudaSpace>]"
${GYSELALIBXX_HOME}/vendor/ddc/include/ddc/chunk.hpp(109): here cannot be referenced -- it is a deleted function
```

Similarly the following error then says that the lambda function itself cannot be created because the function will not work without first capturing the chunk.
```
${GYSELALIBXX_HOME}/vendor/ddc/include/ddc/parallel_for_each.hpp(34): error: function "lambda [](ArgType)->void::<unnamed>(const lambda [](ArgType)->void &)" (declared implicitly) cannot be referenced -- it is a deleted function
```

##  The enclosing parent function for an extended '\_\_host\_\_' '\_\_device\_\_' lambda cannot have private or protected access within its class

Kokkos is bound by restrictions coming from various different low level languages. One of the most surprising of these is the fact that lambdas cannot be found in a private or protected context. In other words we cannot use `ddc::parallel_for_each` or `ddc::parallel_transform_reduce` in a private or protected class method.

In particular a google `TEST` block (as used in all the unit tests) is considered to be a private context. If you wish to put parallel loops in these tests you will need to create a function for your test and call that function from the `TEST` block.

This particular restriction comes from [Cuda](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#extended-lambda-restrictions):

> The enclosing function for the extended lambda must be named and its address can be taken. If the enclosing function is a class member, then the following conditions must be satisfied:
>
> -  All classes enclosing the member function must have a name.
>
> -  The member function must not have private or protected access within its parent class.
>
> -  All enclosing classes must not have private or protected access within their respective parent classes.

## The enclosing parent function for an extended '\_\_host\_\_' '\_\_device\_\_' lambda must allow its address to be taken

Kokkos is bound by restrictions coming from various different low level languages. One of the most surprising of these is the fact that lambdas cannot be found in a function whose address cannot be taken. This includes private and protected class methods (which lead to the previous error message), but also class constructors. If you wish to put parallel loops in a class constructor, you will need to create a function for your loop and call that function from the constructor.

### Error Message

When this rule is not followed you will see the following error message:
```
<PATH_TO_FILE>(<LINE_NUMBER>): error: The enclosing parent function ("<FUNCTION>") for an extended __host__ __device__ lambda cannot have private or protected access within its class
          detected during instantiation of "<CLASS>::<CLASS>(Argtypes)"
```

## A nonstatic member reference must be relative to a specific object

A static function in a class is a function which can be called without having an instance of a class. This can be useful to write functions which return information which is known at compile time (e.g. the degree of a spline) or to write factory functions. A factory function is a function which builds an instance of the class. This is very similar to a constructor but is a little more general. As a result it is capable of returning objects of different types (e.g. subclasses) and running extra calculations before calling the constructor.

For example static functions are sometimes used to initialise classes from a PDI input.

As the static function can be called without having an instance of a class, this error arises when you try to use the members of the class. If the function needs to use the members of the class and not simply call the constructor then the function should not be static.

## X is not defined

This usually occurs when the file containing the definition of `X` is not correctly included. Double check if the `#include <X.hpp>` line is in your file.
