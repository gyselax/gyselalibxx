* We use C++17

# C++ Features
* We don't use raw pointers
* We don't use plain C arrays
* we try not to use macros

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
* type names are CamelCase
* constant & macros use CAPS_NAME
* function and variables use snake_case
* we don't use single letter variables
* we don't rely on case to distinguish between variables
* we name variables with expressive English names (e.g. temperature)
* we don't use variable name conventions from math/physics (u)

# Style
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
