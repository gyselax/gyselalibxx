# C++ Features
* We use C++17
* We don't use raw pointers
* We don't use plain C arrays
* we try not to use macros
* The `auto` keyword should be used for types as `auto&&` and only in one of the following cases:
  - in iterations either in range based for or to refer to iterators
  - to store objects that can not be typed otherwise (e.g. lambdas)
  - to store the result of an expression specifically specifying the type of the generated value

# Parameter passing
* We take out/inout-parameters (those we modify) first
* If there is a single out/inout-parameter, we return it
* For in-parameters (those we use but don't modify)
  - if it's a scalar native type (int, double, ...) we take it by copy
  - if it's a const view type (`std::span<const T>`, `BlockView<T>`, ...) we take it by copy
  - otherwise, we take a cont-ref: `Type const&`
* for out/inout-parameters
  - if it's a modifiable view type (`std::span<T>`, `BlockSpan<T>`, ...) we take it by copy
  - otherwise, we take a ref: `Type&`

# Naming
* we name everything using expressive English names (e.g. `temperature`) and don't use variable
  names from the equations (e.g. `u`, `u_bar_star`)
* files, functions and variables use `snake_case`
* types use `CamelCase`
* macros use `ALL_CAPS`
* non-static member variables names begin with an `m_` prefix
* static member variables names begin with an `s_` prefix
* we don't use single letter variables
* we don't rely on case to distinguish between variables

# Style
* We use the style specified by the `.clang-format` file using clang-format 10
* we do not use numerical values in the code except to initialize a named constexpr documenting
  the semantic of the value
* In a class
  - we put all member types first (public, then protected, then private),
  - followed by static member variables (public, then protected, then private),
  - followed by non-static member variables (public, then protected, then private),
  - followed by static member functions (public, then protected, then private),
  - followed by non-static member functions (public, then protected, then private),
    * the constructors first
    * then the destructor
    * then the various operators
    * then the accessors
    * then the more complex functions
* We comment our code with Doxygen
* We use at @keywords in Doxygen
* we use east-const: `int const` rather than `const int`

# Operators

## Interfaces

We define an *interface* to be an empty class (i.e. without any member
variable) offering a pure virtual call operator: `operator()` and can
offer various overloads of this operator with different parameters but
similar behaviour.

Interfaces should:
* be prefixed by `I` as in `IVlasovSolver`,
* explicitly define a virtual destructor,
* implicitly define constructors and assignment operators.

## Implementation

Implementation classes should implement at least one interface.
When relying on other operations, these classes should use the dependency
injection pattern.
Each implementation class should take `const` references to the used operations
interfaces as parameters of the constructor and store them so as to call them in
the `operator()` implementation.
