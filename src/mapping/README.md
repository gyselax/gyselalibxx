# Mappings

This folder contains code describing tools for handling different coordinate systems.

The current coordinate transformations implemented are:

- Analytically invertible coordinate transformations such as:
	-  Circular coordinate transformation (CircularToCartesian/CartesianToCircular):
		-  $x(r,\theta) = r \cos(\theta),$
	 	-  $y(r,\theta) = r \sin(\theta).$
	-  Czarny coordinate transformation (CzarnyToCartesian/CartesianToCzarny):
		-  $x(r,\theta) = \frac{1}{\epsilon} \left( 1 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \right),$
		-  $y(r,\theta) = \frac{e\xi r \sin(\theta)}{2 -\sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} },$
		 with $\xi = 1/\sqrt{1 - \epsilon^2 /4}$ and $e$ and $\epsilon$ given as parameters.
    - Barycentric coordinate transformation (BarycentricToCartesian):
        -  $`(c_1, c_2, c_3) -> (x,y)`$
- Discrete coordinate transformation defined on B-splines (DiscreteToCartesian):
	-  $`x(r,\theta) = \sum_k c_{x,k} B_k(r,\theta),`$
	-  $`y(r,\theta) = \sum_k c_{y,k} B_k(r,\theta).`$
- Combined coordinate transformation which combines two of the coordinate transformations above.

The tools are:
- InverseJacobianMatrix : this tool calculates the inverse Jacobian matrix on the specified coordinate system.
- InvJacobianOPoint : this tool calculates the inverse Jacobian matrix at the O-point on the specified coordinate system.
- MetricTensor : this tool calculates the metric tensor associated with a coordinate transformation.
- VectorMapper : this tool helps when converting vectors stored in a `VectorField` from one coordinate system to another.
- other static analysis tools found in `mapping_tools.hpp`

