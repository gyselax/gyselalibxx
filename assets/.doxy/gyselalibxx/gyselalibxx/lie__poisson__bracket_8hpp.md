

# File lie\_poisson\_bracket.hpp



[**FileList**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**lie\_poisson\_bracket.hpp**](lie__poisson__bracket_8hpp.md)

[Go to the source code of this file](lie__poisson__bracket_8hpp_source.md)



* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "gradient.hpp"`
* `#include "metric_tensor_evaluator.hpp"`
* `#include "static_tensors.hpp"`
* `#include "vector_field.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**LiePoissonBracket**](classLiePoissonBracket.md) &lt;class Mapping3D, class MappingCoord&gt;<br>_A class which implements a gyrokinetic Poisson bracket operator. The implemented equation is:_ \(\{F, G\} = b\dot(\nabla F \cross \nabla G)\) _with_\(b= \mathbf{B} / B\) _the unitary magnetic field, i.e:_\(\{F, G\} = {\cal J}_{\rm x}^{-1}\epsilon^{ijk}\partial_{x^i} F \partial_{x^j} G b_k\) _with_\({\cal J}_{\rm x}\) _the jacobian of the system,_\(b_k\) _the covariant components of b and_\(\epsilon^{ijk}\) _the Levi-Civita symbol._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/lie_poisson_bracket.hpp`

