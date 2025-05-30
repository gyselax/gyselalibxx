# Coding Standards

## C++ Features

- We use C++17
- We don't use raw pointers
- We don't use plain C arrays
- We try not to use macros
- We try to `static_assert` all assumed information that can be known at compile-time
- We use `std::runtime_error` on CPU and `Kokkos::abort` on GPU to enforce suppositions which are dependent on user input. These checks should appear in constructors
- We use `assert` on CPU and `KOKKOS_ASSERT` on GPU to enforce suppositions which are independent of user input.
- The `auto` keyword should be avoided where possible. It is permitted in one of the following cases:
  - to store objects that can not be typed otherwise (e.g. lambdas)
  - to store the result of an expression specifically specifying the type of the generated value

## Parameter passing

- We take out/inout-parameters (those we modify) first
- For in-parameters (those we use but don't modify)
  - if it's a scalar native type (int, double, ...) we take it by copy
  - if it's a const view type (`ConstField<...>`, `DVectorConstField<...>`, ...) we take it by copy
  - otherwise, we take a const-ref: `Type const&`
- For inout-parameters
  - if it's a modifiable view type (`Field<...>`, `DVectorField<...>`, ...) we take it by copy
  - otherwise, we take a ref: `Type&`

## Naming

- Files, functions and variables use `snake_case`.
- Types use `CamelCase`.
- Macros use `ALL_CAPS`.
- Non-static member variables names begin with an `m_` prefix.
- Static member variables names begin with an `s_` prefix.
- We don't use single letter variables.
- We don't rely on case to distinguish between variables.
- There are two types of DDC objects representing a multidimensional array : `FieldMem` (which possesses the data) and `Field` (which does not own the data but can be captured by `KOKKOS_LAMBDA`). We suffix `FieldMem` objects with `_alloc` if both variables are needed locally.
- If a variable is mirrored between host (CPU) and device (GPU) memories, the variable representing data on host is `_host` suffixed.
- Capturing classes members through `KOKKOS_LAMBDA` or `KOKKOS_CLASS_LAMBDA` may be complicated, we often need to copy-by-reference the member to a local variable, which must be `_proxy` suffixed.
- Types representing continuous dimensions ($x$, $y$, $r$, $\theta$, etc) are named similarly to the mathematical dimension (`X`, `Y`, `R`, `Theta`, etc).
- Types representing a discretised grid along a continuous dimension are named `GridDim` where `Dim` is the continuous dimension.
- Types representing B-splines along a continuous dimension are named `BSplinesDim` where `Dim` is the continuous dimension.
- Grids on an unknown dimension are named `Grid1D`.
- Aliases for `Idx<..>`, `IdxStep<..>`, `IdxRange<..>` should begin with the name of the template (e.g. `IdxX`).

## Style

- We use the style specified by the `.clang-format` file using clang-format 10.
- we do not use numerical values in the code except to initialise a named constexpr documenting
  the semantic of the value.
- In a class
  - we put all member types first (public, then protected, then private),
  - followed by static member variables (public, then protected, then private),
  - followed by non-static member variables (public, then protected, then private),
  - followed by static member functions (public, then protected, then private),
  - followed by non-static member functions (public, then protected, then private),
    - the constructors first
    - then the destructor
    - then the various operators
    - then the accessors
    - then the more complex functions
- We comment our code with Doxygen.
- We use at @@keywords in Doxygen.
- we use east-const: `int const` rather than `const int`.

## GPU management

- GPU code is the default.
- When calling methods beginning with the prefix `parallel_` we always specify the execution space.
- When calling mirror functions, the host execution space is never passed as an argument.

## Operators

### Interfaces

We define an *interface* to be an empty class (i.e. without any member
variable) offering a pure virtual call operator: `operator()` and can
offer various overloads of this operator with different parameters but
similar behaviour.

Interfaces should:

- be prefixed by `I` as in `IVlasovSolver`,
- explicitly define a virtual destructor,
- implicitly define constructors and assignment operators.

### Implementation

Implementation classes should implement at least one interface.
When relying on other operations, these classes should use the dependency
injection pattern.
Each implementation class should take `const` references to the used operations
interfaces as parameters of the constructor and store them so as to call them in
the `operator()` implementation.

## Code Organisation

### Object-Oriented vs. Functional Programming

There are multiple coding paradigms possible when writing code in C++. Two of the most common are object-oriented programming and functional programming. In object-oriented programming the aim is to make code more readable by grouping concepts within a class. Therefore functions are saved with the data on which they operate. On the other hand, in functional programming the aim is to make code more readable by making it more obvious where data is modified. This is done by separating data storage from operators. This way, data can never be modified inside its class without passing through an operator. In our code the functional programming paradigm was chosen as it allows us to write the code in a way which more closely resembles the equations.

The code is therefore generally split into two types of classes: Operators and Data Storage.
In order to ensure that you are following the functional programming paradigm when writing your code you should try to respect this separation.

#### Operators

An operator is an object which can be described by an equation. For example a Poisson solver, an advection operator or a spline interpolation.

Classes describing operators should contain as little data as possible. Data should only be stored in the operator if it is not relevant anywhere else in the code. If removing this operator to replace it with a new method would result in this data also being removed then it may be relevant to store it in the class.

All operators should implement the function `operator()`. Usually this is the only function in the operator other than the constructor.

If the operator relies on other operators (e.g. a semi-Lagrangian advection which requires an interpolation operator) then this dependency should be passed to the operator's constructor. It should be stored as a reference inside the function.

#### Data Storage

It is rare to need to actually write a data storage class. In most cases they can be defined with a "using-declaration" using DDC types. This ensures that the data is defined on the relevant dimensions which in turn reduces the chance of errors arising from operating on the wrong dimensions of an array. For more information see [DDC](https://github.com/CExA-project/ddc).

### Template Strategy

Templating is a very useful tool to improve performance, however it increases compilation overhead and reduces code readability, especially for developers less familiar with C++. We therefore try to strike a balance with the use of templates. Templates are therefore recommended in the following cases:

- In performance bottlenecks to improve results.
- Where they significantly decrease code duplication (~2 additional copies avoided).
