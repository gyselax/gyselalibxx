# Interface derivatives test

## Single interface derivatives calculator

Test `SingleInterfaceDerivativesCalculator` with

- Hermite boundary conditions and additional interpolation points as closure condition.
- Uniform and non-uniform meshes.
- Different patch connections as follows:

![Illustration test 1](../../../docs/images/interface_derivatives/fig\_test\_1.png "")

![Illustration test 2](../../../docs/images/interface_derivatives/fig\_test\_2.png "")

![Illustration test 3](../../../docs/images/interface_derivatives/fig\_test\_3.png "")

![Illustration test 4](../../../docs/images/interface_derivatives/fig\_test\_4.png "")

![Illustration test 5](../../../docs/images/interface_derivatives/fig\_test\_5.png "")

with $\theta$ and $\xi$ periodic.

## Interface derivative matrix with exact formula

Test `InterfaceExactDerivativeMatrix` with the following test cases

- In `interface_derivative_matrix_Greville_periodic_test`: we test with additional interpolation points as closure condition for the $y$-axis(`ddc::BoundCond::GREVILLE`) at the North and the South of the global domain, and periodic boundary conditions for the $x$-axis (`ddc::BoundCond::PERIODIC`) at the West and the East of the global domain. All the patches follow the same orientation as the global domain.

![Illustration test Greville and periodic boundary conditions](../../../docs/images/interface_derivatives/fig5\_example\_9\_patches.png "")

- Hermite boundary conditions with the following layouts to test the signs of the derivatives.

  - In `interface_derivative_matrix_REVERSE_PATCH1_test`: we test that the correct values and signs are selected, especially for the left boundary.
![Illustration test with reversed Patch1](../../../docs/images/interface_derivatives/fig7\_test\_REVERSE\_PATCH1.png "")

  - In `interface_derivative_matrix_REVERSE_PATCH2_test`: we test that the correct values and signs are selected, especially after several inversions.
![Illustration test with reversed Patch2](../../../docs/images/interface_derivatives/fig8\_test\_REVERSE\_PATCH2.png "")

  - In `interface_derivative_matrix_REVERSE_PATCH3_test`: we test that the correct values and signs are selected, especially for the right boundary.
![Illustration test with reversed Patch3](../../../docs/images/interface_derivatives/fig9\_test\_REVERSE\_PATCH3.png "")

  - In `interface_derivative_matrix_CHANGE_BOUND1_test`: we that for an agreement of the direction of the derivatives along $`x_1`$ but not on $`y_1`$
    for the left boundary.
![Illustration test with changed bound Patch1](../../../docs/images/interface_derivatives/fig10\_test\_CHANGE\_BOUND1.png "")

  - In `interface_derivative_matrix_CHANGE_BOUND3_test`: we that for an agreement of the direction of the derivatives along $`y_3`$ but not on $`x_3`$
    for the right boundary.
![Illustration test with changed bound Patch3](../../../docs/images/interface_derivatives/fig11\_test\_CHANGE\_BOUND3.png "")
