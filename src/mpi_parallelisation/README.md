# Parallelisation

This folder contains all code related to the MPI parallelisation. The parallelisation is mainly managed by two concepts:
- Layout
- TransposeOperator

This folder contains the implementation of the superclasses for these two concepts as well as some specific implementations for these concepts. There are also some tools which are useful for the implementation of these concepts.

## Layout

The layouts describe how the data should be laid out in memory across different MPI processes. It indicates the normal ordering of dimensions in the domain associated with the layout and which of those dimensions may contains data distributed across multiple MPI processes.
These classes implement an operator which allows the local domain to be obtained from a global domain and specifications about the MPI configuration (number of ranks and current rank).

## TransposeOperator

The transpose operators are the operators which are used to move from one layout to another. They send and receive data between MPI processes.

## Alltoall Transpose Operator

The alltoall transpose operator is based on the transpose operator present in the Fortran version of Gysela. It uses MPI's Alltoall operator to move from a layout distributed over a given set of dimensions to another layout distributed over an orthogonal set of dimensions. This is achieved by reordering the data such that the data blocks to be sent to each MPI rank are contiguous. Finally after the Alltoall call the data is reordered back into the expected final layout.

### Example
Let us consider the 5D domain (Sp, R, Theta, Vpar, Mu) with the following number of points in each dimension : $`(n_{sp} = 2, n_r = 4, n_\theta = 16, n_{vpar} = 8, n_\mu = 4)`$.

Suppose that we distribute this data in two different ways:
1.  **PoloidalLayout** - Distributed in $(Vpar, Mu)$.
2.  **CollisionalLayout** - Distributed in $(R, Theta)$.

If we have 16 MPI ranks then each rank has a domain with the following number of points in each dimension for each layout:
1.  **PoloidalLayout** - $`(n_{sp} = 2, n_{vpar} = 1, n_\mu = 2, n_r = 4, n_\theta = 16)`$
2.  **CollisionalLayout** - $`(n_{sp} = 2, n_r = 1, n_\theta = 4, n_{vpar} = 8, n_\mu = 4)`$

Suppose that the alltoall transpose operator is called with an input type of `DField<IdxRange<Species, GridVpar, GridMu, GridR, GridTheta>>` on the poloidal layout and an output type of `DField<IdxRange<Species, GridR, GridTheta, GridVpar, GridMu>>` on the collisional layout.
The steps carried out by the alltoall transpose operator are:

1.  Break the field to introduce the splits along which the field will be split. This operation does not require any copying.
    - Before:
        - **type**: `DField<IdxRange<Species, GridVpar, GridMu, GridR, GridTheta>>`
        - **size**: (2, 1, 2, 4, 16)
    - After:
        - **type**: `DField<IdxRange<Species, GridVpar, GridMu, MPIDim<GridR>, GridR, MPIDim<GridTheta>, GridTheta>>`
        - **size**: (2, 1, 2, 4, 1, 4, 4)

2.  Reorder the dimensions so the MPI dimensions are grouped in the first dimensions of the domain:
    - Before:
        - **type**: `DField<IdxRange<Species, GridVpar, GridMu, MPIDim<GridR>, GridR, MPIDim<GridTheta>, GridTheta>>`
        - **size**: (2, 1, 2, 4, 1, 4, 4)
    - After:
        - **type**: `DField<IdxRange<MPIDim<GridR>, MPIDim<GridTheta>, Species, GridVpar, GridMu, GridR, GridTheta>>`
        - **size**: (4, 4, 2, 1, 2, 1, 4)

3.  Call AlltoAll to send the subblocks `field[mpi_idx_r, mpi_idx_theta]` to each MPI rank. The data received from each rank will allow the Vpar and Mu dimensions to be reassembled.
    - Before:
        - **type**: `DField<IdxRange<MPIDim<GridR>, MPIDim<GridTheta>, Species, GridVpar, GridMu, GridR, GridTheta>>`
        - **size**: (4, 4, 2, 1, 2, 1, 4)
    - After:
        - **type**: `DField<IdxRange<MPIDim<GridVpar>, MPIDim<GridMu>, Species, GridVpar, GridMu, GridR, GridTheta>>`
        - **size**: (4, 4, 2, 1, 2, 1, 4)

4.  Reorder the dimensions so data relative to the Vpar and Mu dimensions are placed contiguously and dimensions are ordered as expected in the output layout:
    - Before:
        - **type**: `DField<IdxRange<MPIDim<GridVpar>, MPIDim<GridMu>, Species, GridVpar, GridMu, GridR, GridTheta>>`
        - **size**: (4, 4, 2, 1, 2, 1, 4)
    - After:
        - **type**: `DField<IdxRange<Species, GridR, GridTheta, MPIDim<GridVpar>, GridVpar, MPIDim<GridMu>, GridMu>>`
        - **size**: (2, 1, 4, 4, 2, 4, 1)

5.  Remove the metadata describing the MPI splits. This operation does not require any copying.
    - Before:
        - **type**: `DField<IdxRange<Species, GridR, GridTheta, MPIDim<GridVpar>, GridVpar, MPIDim<GridMu>, GridMu>>`
        - **size**: (2, 1, 4, 4, 2, 4, 1)
    - After:
        - **type**: `DField<IdxRange<Species, GridR, GridTheta, GridVpar, GridMu>>`
        - **size**: (2, 1, 4, 8, 4)
