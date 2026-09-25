# Interface derivatives test

## Interface derivative coefficients

Test `InterfaceDerivCoeffs` with

- Hermite boundary conditions and additional interpolation points as closure condition.
- Uniform and non-uniform meshes.
- Different patch connections as follows:

![Illustration test 1](../../../docs/images/interface_derivatives/fig\_test\_1.png "")

![Illustration test 2](../../../docs/images/interface_derivatives/fig\_test\_2.png "")

![Illustration test 3](../../../docs/images/interface_derivatives/fig\_test\_3.png "")

![Illustration test 4](../../../docs/images/interface_derivatives/fig\_test\_4.png "")

![Illustration test 5](../../../docs/images/interface_derivatives/fig\_test\_5.png "")

with $\theta$ and $\xi$ periodic.

## Interface derivative matrix with approximation formula

Test `InterfacesDerivativeCalculator` with the following test cases

- In `interfaces_derivative_matrix_Greville_periodic_test`: we test with additional interpolation points as closure condition for the $y$-axis(`ddc::SplineBuilderClosure::GREVILLE`) at the North and the South of the global domain, and periodic boundary conditions for the $x$-axis (`ddc::SplineBuilderClosure::PERIODIC`) at the West and the East of the global domain. All the patches follow the same orientation as the global domain.

![Illustration test Greville and periodic boundary conditions](../../../docs/images/interface_derivatives/fig5\_example\_9\_patches.png "")

- In `interfaces_derivative_matrix_Hermite_test`: we test with Hermite boundary conditions for the following layouts to test the signs of the derivatives.

  - We test that the correct values and signs are selected, especially for the left boundary.

![Illustration test with reversed Patch1](../../../docs/images/interface_derivatives/fig7\_test\_REVERSE\_PATCH1.png "")

  - We test that the correct values and signs are selected, especially after several inversions.

![Illustration test with reversed Patch2](../../../docs/images/interface_derivatives/fig8\_test\_REVERSE\_PATCH2.png "")

  - We test that the correct values and signs are selected, especially for the right boundary.

![Illustration test with reversed Patch3](../../../docs/images/interface_derivatives/fig9\_test\_REVERSE\_PATCH3.png "")

  - We test with an agreement of the direction of the derivatives along $`x_1`$ but not on $`y_1`$
    for the left boundary.

![Illustration test with changed bound Patch1](../../../docs/images/interface_derivatives/fig10\_test\_CHANGE\_BOUND1.png "")

  - We test with an agreement of the direction of the derivatives along $`y_3`$ but not on $`x_3`$
    for the right boundary.

![Illustration test with changed bound Patch3](../../../docs/images/interface_derivatives/fig11\_test\_CHANGE\_BOUND3.png "")
