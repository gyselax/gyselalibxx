

# Struct PolarBSplinesRTheta



[**ClassList**](annotated.md) **>** [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md)








Inherits the following classes: [PolarBSplines](classPolarBSplines.md)
















## Public Types inherited from PolarBSplines

See [PolarBSplines](classPolarBSplines.md)

| Type | Name |
| ---: | :--- |
| typedef [**BSplinesR**](structBSplinesR.md) | [**BSplinesR\_tag**](classPolarBSplines.md#typedef-bsplinesr_tag)  <br>_The radial bspline from which the polar B-splines are constructed._  |
| typedef [**BSplinesTheta**](structBSplinesTheta.md) | [**BSplinesTheta\_tag**](classPolarBSplines.md#typedef-bsplinestheta_tag)  <br>_The poloidal bspline from which the polar B-splines are constructed._  |
| typedef typename BSplinesR::continuous\_dimension\_type | [**R**](classPolarBSplines.md#typedef-r)  <br>_The tag for the radial direction of the B-splines._  |
| typedef typename BSplinesTheta::continuous\_dimension\_type | [**Theta**](classPolarBSplines.md#typedef-theta)  <br>_The tag for the poloidal direction of the B-splines._  |
| typedef [**PolarBSplines**](classPolarBSplines.md) | [**discrete\_dimension\_type**](classPolarBSplines.md#typedef-discrete_dimension_type)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_idx\_range\_type**](classPolarBSplines.md#typedef-tensor_product_idx_range_type)  <br> |
| typedef IdxStep&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_idx\_step\_type**](classPolarBSplines.md#typedef-tensor_product_idx_step_type)  <br> |
| typedef Idx&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type)  <br> |












## Public Static Attributes inherited from PolarBSplines

See [PolarBSplines](classPolarBSplines.md)

| Type | Name |
| ---: | :--- |
|  int constexpr | [**continuity**](classPolarBSplines.md#variable-continuity)   = `C`<br>_The continuity enforced by the B-splines at the singular point._  |
































## Public Static Functions inherited from PolarBSplines

See [PolarBSplines](classPolarBSplines.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**get\_2d\_index**](classPolarBSplines.md#function-get_2d_index) (Idx&lt; DDim &gt; const & idx) <br> |
|  KOKKOS\_FUNCTION Idx&lt; DDim &gt; | [**get\_polar\_index**](classPolarBSplines.md#function-get_polar_index) ([**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) const & idx) <br> |
|  KOKKOS\_FUNCTION Field&lt; ElementType, IdxRange&lt; DDim &gt;, MemorySpace &gt; | [**get\_singular\_subset**](classPolarBSplines.md#function-get_singular_subset) (Field&lt; ElementType, IdxRange&lt; DDim &gt;, MemorySpace &gt; coeffs) <br>_Get the subset of the coefficients associated with singular basis splines._  |
|  KOKKOS\_FUNCTION Field&lt; ElementType, [**tensor\_product\_idx\_range\_type**](classPolarBSplines.md#typedef-tensor_product_idx_range_type), MemorySpace &gt; | [**get\_tensor\_product\_subset**](classPolarBSplines.md#function-get_tensor_product_subset) (Field&lt; ElementType, IdxRange&lt; DDim &gt;, MemorySpace &gt; coeffs) <br>_Get the subset of the coefficients of tensor product basis splines._  |
|  constexpr std::size\_t | [**n\_singular\_basis**](classPolarBSplines.md#function-n_singular_basis) () <br> |
|  constexpr KOKKOS\_FUNCTION IdxRange&lt; DDim &gt; | [**singular\_idx\_range**](classPolarBSplines.md#function-singular_idx_range) () <br>_Get the IdxRange containing the indices of the b-splines which traverse the singular point._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

