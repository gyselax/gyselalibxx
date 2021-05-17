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
  - if it's a view type (`std::span<const T>`, `BlockView<const T>`, ...) we take it by copy
  - otherwise, we take a cont-ref: `Type const&`
* for out/inout-parameters
  - if it's a view type (`std::span<T>`, `BlockView<T>`, ...) we take it by copy
  - otherwise, we take a ref: `Type&`

# Naming
* type names *must* use CamelCase
* constant & macro names *must* be ALL_CAPS
* function and variables names *must* use snake_case
* non-static member variables names *should* begin with an `m_` prefix
* static member variables names *should* begin with an `s_` prefix
* we don't use single letter variables
* we don't rely on case to distinguish between variables
* we name variables with expressive English names (e.g. temperature)
* we don't use variable name conventions from math/physics (u)

# Style
* We *must* apply the style specified by the `.clang-format` file using clang-format 10
* we *must not* use numerical values in the code except to initialize a named constexpr documenting
  the semantic of the value
* In a class
  - we put all member types first (public, then protected, then private),
  - then all member variables (public, then protected, then private),
  - then all static member functions (public, then protected, then private),
  - then all non-static member functions (public, then protected, then private),
    * the constructors first
    * then the destructor
    * then the various operators
    * then the accessors
    * then the more complex functions
* We comment our code with Doxygen
* We use at @keywords in Doxygen
* we use east-const (yoda-const): `int const` rather than `const int`
