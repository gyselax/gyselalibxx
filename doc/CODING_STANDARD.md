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

# Style
* we use post-const (yoda-const): `int const` rather than `const int`
