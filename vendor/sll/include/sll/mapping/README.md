# Mappings

The `mapping/` folder containts the code defining the mappings used in simulations.
The current mappings implemented are:

- Analytically invertible mappings such as:
	-  Circular mapping (CircularToCartesian):
		-  $x(r,\theta) = r \cos(\theta),$
	 	-  $y(r,\theta) = r \sin(\theta).$
	-  Czarny mapping (CzarnyToCartesian):
		-  $x(r,\theta) = \frac{1}{\epsilon} \left( 1 - \sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} \right),$
		-  $y(r,\theta) = \frac{e\xi r \sin(\theta)}{2 -\sqrt{1 + \epsilon(\epsilon + 2 r \cos(\theta)} },$
		 with $\xi = 1/\sqrt{1 - \epsilon^2 /4}$ and $e$ and $\epsilon$ given as parameters.
- Discrete mappings defined on bsplines (DiscreteToCartesian):
	-  $`x(r,\theta) = \sum_k c_{x,k} B_k(r,\theta),`$
	-  $`y(r,\theta) = \sum_k c_{y,k} B_k(r,\theta).`$
- Refined discrete mappings defined on bsplines (RefinedDiscreteToCartesian): a version of DiscreteToCartesian mapping built on a refined grid.  
