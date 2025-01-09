# Matrix tools

This folder describes various matrix classes which can be used to store matrices in order to solve matrix equations.

They can be approximately grouped into 2 groups:
1. Classes inheriting from `Matrix`. These classes solve matrix equations using LAPACK on CPU. These matrices should be created using the factory methods provided in the `Matrix` class.
2. Classes inheriting from `MatrixBatch`. These classes can solve matrix equations on GPU. They solve 1D equations batched over another dimension.
