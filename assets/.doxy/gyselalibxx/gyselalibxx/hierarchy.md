
# Class Hierarchy

This inheritance list is sorted roughly, but not completely, alphabetically:


* **class** [**AdvectionFieldFinder**](classAdvectionFieldFinder.md) _Solve the Poisson-like equation and return the electric field for the coupled Vlasov equation._ 
* **class** [**BarycentricToCartesian**](classBarycentricToCartesian.md) _A class to convert barycentric coordinates to Cartesian coordinates on a triangle._ 
* **class** [**BslAdvection1D**](classBslAdvection1D.md) _A class which computes the advection along the dimension of interest GridInterest._ 
* **class** [**BslAdvectionPolar**](classBslAdvectionPolar.md) _Define an advection operator on 2D_  _domain._
* **class** [**IAdvectionSpatial**](classIAdvectionSpatial.md) _A class which provides an advection operator._     
    * **class** [**BslAdvectionSpatial**](classBslAdvectionSpatial.md) _A class which computes the spatial advection along the dimension of interest_ [_**GridX**_](structGridX.md) _. Working for every Cartesian geometry._
* **class** [**IAdvectionVelocity**](classIAdvectionVelocity.md) _A class which provides an advection operator._     
    * **class** [**BslAdvectionVelocity**](classBslAdvectionVelocity.md) _A class which computes the velocity advection along the dimension of interest GridV. Working for every Cartesian geometry._ 
* **class** [**ITimeSolverRTheta**](classITimeSolverRTheta.md) _Base class for the time solvers._     
    * **class** [**BslExplicitPredCorrRTheta**](classBslExplicitPredCorrRTheta.md) _A second order explicit predictor-corrector for the Vlasov-Poisson equations._ 
    * **class** [**BslImplicitPredCorrRTheta**](classBslImplicitPredCorrRTheta.md) _A second order implicit predictor-corrector for the Vlasov-Poisson equations._ 
    * **class** [**BslPredCorrRTheta**](classBslPredCorrRTheta.md) _Predictor-corrector for the Vlasov-Poisson equations._ 
* **class** [**IEquilibrium**](classIEquilibrium.md) _An abstract class for initialising a distribution function in (species,vpar,mu)._     
    * **class** [**BumpontailEquilibrium**](classBumpontailEquilibrium.md) _A class that initialises the distribution function as a sum of two Maxwellian functions._ 
    * **class** [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) _Equilibrium operator as Maxwellian. This initialises all species._ 
    * **class** [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) _Equilibrium operator as Maxwellian. This initialises all species._ 
    * **class** [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) _Equilibrium operator as Maxwellian. This initialises all species._ 
* **class** [**CartesianToBarycentric**](classCartesianToBarycentric.md) _A class to convert Cartesian coordinates to barycentric coordinates on a triangle._ 
* **class** [**CartesianToCircular**](classCartesianToCircular.md) _A class for describing the circular 2D mapping._ 
* **class** [**CartesianToCylindrical**](classCartesianToCylindrical.md) _A class for describing the cylindrical 3D mapping._ 
* **class** [**CartesianToCzarny**](classCartesianToCzarny.md) _A class for describing the Czarny 2D mapping._ 
* **class** [**IPartialDerivative**](classIPartialDerivative.md) _An abstract class for a partial derivative operator._     
    * **class** [**CentralFDMPartialDerivative**](classCentralFDMPartialDerivative.md) _A class which implements a partial derivative operator using a finite differences calculation of order two. A decentered scheme is used at the boundary, whereas centred finite difference are used inside the domain._ 
* **class** [**IPartialDerivativeCreator**](classIPartialDerivativeCreator.md) _An abstract class which provides a create\_instance function to instantiate an object of the_ [_**IPartialDerivative**_](classIPartialDerivative.md) _class where required._    
    * **class** [**CentralFDMPartialDerivativeCreator**](classCentralFDMPartialDerivativeCreator.md) _A class which stores information necessary to create a pointer to an instance of the_ [_**CentralFDMPartialDerivative**_](classCentralFDMPartialDerivative.md) _class._
* **class** [**IChargeDensityCalculator**](classIChargeDensityCalculator.md) _A class which computes charges density._     
    * **class** [**ChargeDensityCalculator**](classChargeDensityCalculator.md) _A class which computes charges density with Kokkos._ 
    * **class** [**ChargeDensityCalculator**](classChargeDensityCalculator.md) _A class which computes charges density with Kokkos._ 
    * **class** [**MpiChargeDensityCalculator**](classMpiChargeDensityCalculator.md) _A class which computes charges density with Kokkos._ 
* **class** [**CircularToCartesian**](classCircularToCartesian.md) _A class for describing the circular 2D mapping._ 
* **class** [**CollisionConfiguration**](classCollisionConfiguration.md) _Class to collect information to initialise the collision operator for a SpVparMu geometry._ 
* **class** [**CollisionOperator**](classCollisionOperator.md) _A class which computes the collision operator in (Sp,vpar,mu)._ 
* **class** [**IRightHandSide**](classIRightHandSide.md) _An abstract class representing a source in Boltzmann equation._     
    * **class** [**CollisionsInter**](classCollisionsInter.md) _Class describing the inter-species collision operator._ 
    * **class** [**CollisionsIntra**](classCollisionsIntra.md) _Class describing the intra-species collision operator._ 
    * **class** [**KineticSource**](classKineticSource.md) _A class that describes a source of particles._ 
    * **class** [**KrookSourceAdaptive**](classKrookSourceAdaptive.md) _A class that describes a source of particles._ 
    * **class** [**KrookSourceConstant**](classKrookSourceConstant.md) _A class that describes a source of particles._ 
* **class** [**CombinedMapping**](classCombinedMapping.md) _A class which describes a mapping which is constructed by combining two mappings. Let us denote Mapping1 as_  _and Mapping2 as_ _then this mapping represents:_ _._
* **class** [**ITimeStepper**](classITimeStepper.md) _The superclass from which all timestepping methods inherit._     
    * **class** [**CrankNicolson**](classCrankNicolson.md) _A class which provides an implementation of a Crank-Nicolson method._ 
    * **class** [**Euler**](classEuler.md) _A class which provides an implementation of an explicit_ [_**Euler**_](classEuler.md) _method._
    * **class** [**RK2**](classRK2.md) _A class which provides an implementation of a second-order Runge-Kutta method._ 
    * **class** [**RK3**](classRK3.md) _A class which provides an implementation of a third-order Runge-Kutta method._ 
    * **class** [**RK4**](classRK4.md) _A class which provides an implementation of a fourth-order Runge-Kutta method._ 
* **class** [**CylindricalToCartesian**](classCylindricalToCartesian.md) _A class for describing the cylindrical 3D mapping._ 
* **class** [**CzarnyToCartesian**](classCzarnyToCartesian.md) _A class for describing the Czarny 2D mapping._ 
* **class** [**DerivField**](classDerivField.md) 
* **class** [**DerivFieldCommon**](classDerivFieldCommon.md) 
* **class** [**DerivFieldCommon&lt; FieldType, IdxRange&lt; DDims... &gt; &gt;**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md) _An abstract class which holds a chunk of memory describing a field and its derivatives. This is the superclass for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _and_[_**DerivField**_](classDerivField.md) _._
* **class** [**DerivFieldMem**](classDerivFieldMem.md) 
* **class** [**DiocotronDensitySolution**](classDiocotronDensitySolution.md) _The diocotron exact solution of the density_  _._
* **class** [**DiscreteToCartesian**](classDiscreteToCartesian.md) _A class for describing discrete 2D mappings from the logical domain to the physical domain._ 
* **class** [**DiscreteToCartesianBuilder**](classDiscreteToCartesianBuilder.md) _A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._
* **class** [**EdgeTransformation**](classEdgeTransformation.md) _Transform a coordinate or an index from one edge to the one on the other edge._ 
* **class** [**IPoissonSolver**](classIPoissonSolver.md) 
* **class** [**FFTPoissonSolver**](classFFTPoissonSolver.md) 
* **class** [**FluidMoments**](classFluidMoments.md) _A class that computes fluid moments of the distribution function._ 
* **class** [**GaussLegendre**](classGaussLegendre.md) _An operator for constructing a Gauss-Legendre quadrature._ 
* **class** [**GeometryVxVyXY**](classGeometryVxVyXY.md) _A class providing aliases for useful subindex ranges of the geometry when the data is saved with the velocity dimensions distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections._ 
* **class** [**GeometryXVx**](classGeometryXVx.md) _A class providing aliases for useful subindex ranges of the geometry. It is used as template parameter for generic dimensionality-agnostic operators such as advections._ 
* **class** [**GeometryXYVxVy**](classGeometryXYVxVy.md) _A class providing aliases for useful subindex ranges of the geometry when the data is saved with the spatial dimensions distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections._ 
* **class** [**Gradient**](classGradient.md) _A class which implements a gradient operator._ 
* **class** [**IBoltzmannSolver**](classIBoltzmannSolver.md) _An abstract class for solving a Boltzmann equation._     
    * **class** [**SplitRightHandSideSolver**](classSplitRightHandSideSolver.md) _A class that solves a Boltzmann equation using Strang's splitting._ 
    * **class** [**SplitVlasovSolver**](classSplitVlasovSolver.md) _A class that solves a Vlasov equation using Strang's splitting._ 
* **class** [**IInitialisation**](classIInitialisation.md) _An abstract class that allows for initialising a distribution function._     
    * **class** [**NoPerturbInitialisation**](classNoPerturbInitialisation.md) _Initialisation operator with no perturbation, i.e the distribution function equal to the Maxwellian._ 
    * **class** [**RestartInitialisation**](classRestartInitialisation.md) _A class that initialises the distribution function from a previous simulation._ 
    * **class** [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) _A class that initialises the distribution function as a perturbed Maxwellian._ 
    * **class** [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) _A class that initialises the distribution function as a perturbed Maxwellian._ 
* **class** [**IInterpolator**](classIInterpolator.md) _A class which provides an interpolating function._ 
* **class** [**IInterpolator2D**](classIInterpolator2D.md) _A class which provides an interpolating function._     
    * **class** [**IPreallocatableInterpolator2D**](classIPreallocatableInterpolator2D.md) _A class which provides access to an interpolating function which can be preallocated where useful._ 
    * **class** [**IPreallocatableInterpolator2D**](classIPreallocatableInterpolator2D.md) _A class which provides access to an interpolating function which can be preallocated where useful._ 
* **class** [**IMPILayout**](classIMPILayout.md) _A super class describing a way in which data may be laid out across MPI processes._ 
* **class** [**IMPITranspose**](classIMPITranspose.md) _A superclass describing an operator for converting from/to different MPI layouts._     
    * **class** [**MPITransposeAllToAll**](classMPITransposeAllToAll.md) _A class describing an operator for converting from/to different MPI layouts using AlltoAll._ 
    * **class** [**MPITransposeAllToAll**](classMPITransposeAllToAll.md) _A class describing an operator for converting from/to different MPI layouts using AlltoAll._ 
* **class** [**IPoissonSolver&lt; IdxRange&lt; ODims... &gt;, IdxRangeFull, MemorySpace, LayoutSpace &gt;**](classIPoissonSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md) 
* **class** [**IPolarFootFinder**](classIPolarFootFinder.md) _Define a base class for all the time integration methods used to find the foot of a characteristic on a polar domain (a polar domain is a domain defined on the_  _plane)._
* **class** [**IQNSolver**](classIQNSolver.md) _Base class for a Quasi-Neutrality solver._     
    * **class** [**NullQNSolver**](classNullQNSolver.md) _Null operator._ 
    * **class** [**NullQNSolver**](classNullQNSolver.md) _Null operator._ 
    * **class** [**QNSolver**](classQNSolver.md) _An operator which solves the Quasi-Neutrality equation using a fast Fourier transform._ 
    * **class** [**QNSolver**](classQNSolver.md) _An operator which solves the Quasi-Neutrality equation using a fast Fourier transform._ 
* **class** [**ITimeSolver**](classITimeSolver.md) _An abstract class for solving a Boltzmann-Poisson system of equations._     
    * **class** [**PredCorr**](classPredCorr.md) _A class that solves a Boltzmann-Poisson system of equations using a predictor-corrector scheme._ 
    * **class** [**PredCorr**](classPredCorr.md) _A class that solves a Boltzmann-Poisson system of equations using a predictor-corrector scheme._ 
* **class** [**IVlasovSolver**](classIVlasovSolver.md) _An abstract class for solving a Vlasov equation._     
    * **class** [**MpiSplitVlasovSolver**](classMpiSplitVlasovSolver.md) _A class that solves a Vlasov equation using Strang's splitting on an MPI distributed mesh._ 
    * **class** [**SplitVlasovSolver**](classSplitVlasovSolver.md) _A class that solves a Vlasov equation using Strang's splitting._ 
* **class** [**IdxRangeSlice**](classIdxRangeSlice.md) _A class which describes a collection of equally spaced Idxs which form a index range._ 
* **class** [**InvJacobianOPoint**](classInvJacobianOPoint.md) _An operator for calculating the inverse of the Jacobian at an O-point. This class is used in_ [_**CombinedMapping**_](classCombinedMapping.md) _to calculate the inverse of the Jacobian at an O-point when one of the mappings does not allow the evaluation of its Jacobian/inverse Jacobian at the O-point._
* **class** [**InvJacobianOPoint&lt; CombinedMapping&lt; CircularToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01CircularToCartesian_3_01R_00_01Theta_00_01X_00_01be6b75c3c69e2165a260584a5fd55276.md) 
* **class** [**InvJacobianOPoint&lt; CombinedMapping&lt; CzarnyToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01CzarnyToCartesian_3_01R_00_01Theta_00_01X_00_01Y_399a54ae9b96ca8e75637deab2a35d95.md) 
* **class** [**InvJacobianOPoint&lt; CombinedMapping&lt; DiscreteToCartesian&lt; X, Y, SplineEvaluator, R, Theta, MemorySpace &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01DiscreteToCartesian_3_01X_00_01Y_00_01SplineEvalu8b276096a791392bbf161b1f2d34e864.md) 
* **class** [**InverseJacobianMatrix**](classInverseJacobianMatrix.md) 
* **class** [**KelvinHelmholtzInstabilityInitialisation**](classKelvinHelmholtzInstabilityInitialisation.md) _Initialise the allfdistribu function._ 
* **class** [**Lagrange**](classLagrange.md) _A class which implements_ [_**Lagrange**_](classLagrange.md) _polynomials._
* **class** [**LeviCivitaTensor**](classLeviCivitaTensor.md) 
* **class** [**MatchingIdxSlice**](classMatchingIdxSlice.md) _Store the conforming indexes of each patch of a given interface._ 
* **class** [**Matrix**](classMatrix.md) _The super class from which matrix classes should inherit. This class is used to solve matrix equations._     
    * **class** [**Matrix\_Banded**](classMatrix__Banded.md) _A matrix class representing a banded matrix._ 
    * **class** [**Matrix\_Corner\_Block**](classMatrix__Corner__Block.md) _A class representing a matrix with the following block pattern:_     
        * **class** [**Matrix\_Centre\_Block**](classMatrix__Centre__Block.md) _A_ [_**Matrix**_](classMatrix.md) _representing a matrix which has a banded region. This matrix must be able to be described by the following block matrices:_
        * **class** [**Matrix\_Periodic\_Banded**](classMatrix__Periodic__Banded.md) _A class representing a periodic banded matrix. A periodic banded matrix is like a banded matrix but additionally contains non- zero values in the corners. I.e it has the following sparsity pattern:_ 
    * **class** [**Matrix\_Dense**](classMatrix__Dense.md) _A class describing a dense matrix._ 
    * **class** [**Matrix\_PDS\_Tridiag**](classMatrix__PDS__Tridiag.md) _A class representing a real symmetric positive definite matrix._ 
* **class** [**MatrixBatch**](classMatrixBatch.md) [_**MatrixBatch**_](classMatrixBatch.md) _superclass for managing a collection of linear systems. The main assumption is that all matrices have the same size. It is also assumed that each matrix is used to solve one equation._    
    * **class** [**MatrixBatchCsr**](classMatrixBatchCsr.md) [_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. This class uses the CSR storage format which needs three arrays, one stores values, the other column indices. The third array contains the count of non-zero inside the matrix lines.(eg:for a given line index i nn\_per\_row[i]= sum of non-zeros until line i) The class returns these arrays (as Kokkos views) with the get\_csr\_views function, it is then possible to fill them outside the class. The sparsity pattern is the same for all matrices, hence column indices are stored only for one system. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._
    * **class** [**MatrixBatchCsr**](classMatrixBatchCsr.md) [_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. This class uses the CSR storage format which needs three arrays, one stores values, the other column indices. The third array contains the count of non-zero inside the matrix lines.(eg:for a given line index i nn\_per\_row[i]= sum of non-zeros until line i) The class returns these arrays (as Kokkos views) with the get\_csr\_views function, it is then possible to fill them outside the class. The sparsity pattern is the same for all matrices, hence column indices are stored only for one system. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._
    * **class** [**MatrixBatchEll**](classMatrixBatchEll.md) [_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. The sparsity pattern is assumed to be the same for all matrices. ie the non-zero components are located at the same places for all matrices. This class uses the ELL storage format which needs two 1D arrays, one stores values the other column indices. The class returns these arrays (as Kokkos views) with the get\_batch\_idx\_and\_vals function, it is then possible to fill them outside the class. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._
    * **class** [**MatrixBatchTridiag**](classMatrixBatchTridiag.md) _A structure for solving a set of independent tridiagonal systems using a direct method. The parallelism operates on the whole collection by dispatching to threads. Each problem is treated sequentially, by the tridiagonal matrix algorithm (TDMA). This solver is stable for tridiagonal matrices which satisfy one of the following conditions:_ 
* **class** [**MetricTensorEvaluator**](classMetricTensorEvaluator.md) _An operator for calculating the metric tensor._ 
* **class** [**Moments**](classMoments.md) [_**Moments**_](classMoments.md) _discrete dimension to access constant attributes related to fluid moments._    
    * **struct** [**GridMom**](structGridMom.md) 
* **class** [**Moments::Impl**](classMoments_1_1Impl.md) [_**Impl**_](classMoments_1_1Impl.md) _object storing attributes in_`MemorySpace` _._
* **class** [**MultipatchConnectivity**](classMultipatchConnectivity.md) _A helper class which provides functionalities to recognise how different patches are connected._ 
* **class** [**MultipatchType**](classMultipatchType.md) _A class to store several objects that are of a type which is templated by the patch._ 
* **class** [**MultipatchSplineBuilder**](classMultipatchSplineBuilder.md) _A class to call all the builders of all the patches once._ 
* **class** [**MultipatchSplineBuilder2D**](classMultipatchSplineBuilder2D.md) _A class to call all the builders of all the patches once._ 
* **class** [**MultipatchSplineEvaluator2D**](classMultipatchSplineEvaluator2D.md) _A class to evaluate all the splines of all the patches at once._ 
* **class** [**OnionPatchLocator**](classOnionPatchLocator.md) [_**Patch**_](structPatch.md) _locator specialised for "onion" geometry._
* **class** [**OnionPatchLocator&lt; MultipatchType&lt; IdxRangeOnPatch, Patches... &gt;, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace &gt;**](classOnionPatchLocator_3_01MultipatchType_3_01IdxRangeOnPatch_00_01Patches_8_8_8_01_4_00_01Logicff6c45b073183ccdfc0de0e4a415a7fa.md) [_**Patch**_](structPatch.md) _locator specialised for "onion" geometry._
* **class** [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) _Type of right-hand side (rhs) function of the Poisson equation._ 
* **class** [**PolarBSplines**](classPolarBSplines.md) 
* **class** [**PolarBSplines::Impl**](classPolarBSplines_1_1Impl.md) 
* **class** [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) _Define an evaluator on polar B-splines._ 
* **class** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) _Define a polar PDE solver for a Poisson-like equation._ 
* **class** [**PredCorrRK2XY**](classPredCorrRK2XY.md) _Predictor-corrector based on_ [_**RK2**_](classRK2.md) _for the guiding-centre model._
* **class** [**Quadrature**](classQuadrature.md) _A class providing an operator for integrating functions defined on known grid points._ 
* **class** [**RefinedDiscreteToCartesianBuilder**](classRefinedDiscreteToCartesianBuilder.md) _A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates an instance which uses more refined splines than the provided builder and evaluator. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._
* **class** [**SpeciesInformation**](classSpeciesInformation.md) [_**Species**_](structSpecies.md) _discrete dimension to access constant attributes related to species._    
    * **struct** [**Species**](structSpecies.md) 
* **class** [**SpeciesInformation::Impl**](classSpeciesInformation_1_1Impl.md) [_**Impl**_](classSpeciesInformation_1_1Impl.md) _object storing attributes in_`MemorySpace` _._
* **class** [**SplineBuilder2DCache**](classSplineBuilder2DCache.md) _A class that stores spline builder coefficients and recomputes them when required._ 
* **class** [**Tensor**](classTensor.md) _A class representing a_ [_**Tensor**_](classTensor.md) _._
* **class** [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md) _A class which evaluates the triangular Bernstein polynomials._ 
* **class** [**TriangularBernsteinPolynomialBasis::Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) 
* **class** [**VectorFieldCommon**](classVectorFieldCommon.md) 
* **class** [**VectorMapper**](classVectorMapper.md) 
* **class** [**VectorMapper&lt; VectorIndexSet&lt; XIn, YIn &gt;, VectorIndexSet&lt; XOut, YOut &gt;, Mapping, ExecSpace &gt;**](classVectorMapper_3_01VectorIndexSet_3_01XIn_00_01YIn_01_4_00_01VectorIndexSet_3_01XOut_00_01YOu77c12468788509067d2c0ef34f5e389c.md) _A class to map vector fields from one coordinate system to another._ 
* **class** [**VortexMergerDensitySolution**](classVortexMergerDensitySolution.md) _Initial condition for the vortex merger simulation._ 
* **class** [**VortexMergerEquilibria**](classVortexMergerEquilibria.md) _Equilibrium solution of a Vlasov-Poissson equations system in polar coordinates._ 
* **class** [**ddcHelper::NonUniformInterpolationPoints**](classddcHelper_1_1NonUniformInterpolationPoints.md) _Helper class for the initialisation of the mesh of interpolation points._ 
* **class** [**tensor\_tools::IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) _A class to capture the description of a tensor indexed at a specific component. This class should not be explicitly declared in user code. It is the output of a call to the index&lt;...&gt; function and is an input to the tensor\_mul function._ 
* **class** [**tensor\_tools::TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md) _A class describing an index of a tensor. For example for a 2x2 metric tensor on an (x,y) plane the element_  _would have the index TensorIndexElement&lt;TensorIndexSetXY, X, X&gt;._
* **struct** [**ConstPolarSpline**](structConstPolarSpline.md) _A structure containing the two ConstFields necessary to define a constant reference to a spline on a set of polar basis splines._ 
* **struct** [**ConstantExtrapolationRuleOnion**](structConstantExtrapolationRuleOnion.md) _Define constant extrapolation rule for onion shape geometries. Struct useful for the MultipatchSplineEvaluator types._  __
* **struct** [**Edge**](structEdge.md) _Define an edge of a given patch._ 
* **struct** [**FEM1DPoissonSolver::FEMSplineEvaluatorBuilder**](structFEM1DPoissonSolver_1_1FEMSplineEvaluatorBuilder.md) 
* **struct** [**FEM1DPoissonSolver::FEMSplineEvaluatorBuilder&lt; IdxRange&lt; DimX... &gt; &gt;**](structFEM1DPoissonSolver_1_1FEMSplineEvaluatorBuilder_3_01IdxRange_3_01DimX_8_8_8_01_4_01_4.md) 
* **struct** [**FluidMoments::MomentDensity**](structFluidMoments_1_1MomentDensity.md) 
* **struct** [**FluidMoments::MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) 
* **struct** [**FluidMoments::MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) 
* **struct** [**GaussLegendreCoefficients**](structGaussLegendreCoefficients.md) _A structure containing the weights and positions associated with a Gauss-Legendre quadrature using NPoints points._ 
* **struct** [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) _An iterator type for the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _._
* **struct** [**IdxRangeToSlice**](structIdxRangeToSlice.md) _A class to create a_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _type from a TypeSeq._
* **struct** [**Interface**](structInterface.md) _Represent a simple sticking of two edges._ 
* **struct** [**MPIDim**](structMPIDim.md) _An internal tag used to dsecribe an artificial dimension describing the MPI rank where the scattered information will be sent to or where the gathered information will be collected from._ 
* **struct** [**Mu**](structMu.md) _Define non periodic magnetic momentum_  _._
* **struct** [**MultipatchSplineBuilder2D::Build\_BuilderType**](structMultipatchSplineBuilder2D_1_1Build__BuilderType.md) 
* **struct** [**MultipatchSplineBuilder2D::Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;**](structMultipatchSplineBuilder2D_1_1Build__BuilderType_3_01Patch_00_01DConstField_3_01IdxRange_3_388990a8744187d12e0f612652c86727.md) 
* **struct** [**MultipatchSplineBuilder::Build\_BuilderType**](structMultipatchSplineBuilder_1_1Build__BuilderType.md) 
* **struct** [**MultipatchSplineBuilder::Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;**](structMultipatchSplineBuilder_1_1Build__BuilderType_3_01Patch_00_01DConstField_3_01IdxRange_3_0143963206948829b7e1024cedaaccd580.md) 
* **struct** [**MultipatchSplineEvaluator2D::eval\_deriv\_type**](structMultipatchSplineEvaluator2D_1_1eval__deriv__type.md) _Tag to indicate that derivative of the spline should be evaluated._ 
* **struct** [**MultipatchSplineEvaluator2D::eval\_type**](structMultipatchSplineEvaluator2D_1_1eval__type.md) _Tag to indicate that the value of the spline should be evaluated._ 
* **struct** [**NullExtrapolationRule**](structNullExtrapolationRule.md) _Define null extrapolation rule common to all geometries._ 
* **struct** [**OutsideEdge**](structOutsideEdge.md) _Define an edge for the outside index range._ [_**OutsideEdge**_](structOutsideEdge.md) _is a pseudo-edge outside the index range used to define interfaces between patches and the outside index range._
* **struct** [**Patch**](structPatch.md) _Base tag for a patch._ 
* **struct** [**Patch&lt; grid1, grid2, bsplines\_dim1, bsplines\_dim2 &gt;**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md) _Tag for a patch._ 
* **struct** [**PolarBSplines::Impl::Corner1Tag**](structPolarBSplines_1_1Impl_1_1Corner1Tag.md) _The tag for the first corner of the Barycentric coordinates._ 
* **struct** [**PolarBSplines::Impl::Corner2Tag**](structPolarBSplines_1_1Impl_1_1Corner2Tag.md) _The tag for the second corner of the Barycentric coordinates._ 
* **struct** [**PolarBSplines::Impl::Corner3Tag**](structPolarBSplines_1_1Impl_1_1Corner3Tag.md) _The tag for the third corner of the Barycentric coordinates._ 
* **struct** [**PolarBSplines::eval\_deriv\_type**](structPolarBSplines_1_1eval__deriv__type.md) 
* **struct** [**PolarBSplines::eval\_type**](structPolarBSplines_1_1eval__type.md) 
* **struct** [**PolarSpline**](structPolarSpline.md) _A structure containing the two Fields necessary to define a reference to a spline on a set of polar basis splines._ 
* **struct** [**PolarSplineEvaluator::eval\_deriv\_r\_theta\_type**](structPolarSplineEvaluator_1_1eval__deriv__r__theta__type.md) _Tag for the evaluation of the cross derivative of the function._ 
* **struct** [**PolarSplineEvaluator::eval\_deriv\_r\_type**](structPolarSplineEvaluator_1_1eval__deriv__r__type.md) _Tag for the evaluation of the derivative on the first dimension._ 
* **struct** [**PolarSplineEvaluator::eval\_deriv\_theta\_type**](structPolarSplineEvaluator_1_1eval__deriv__theta__type.md) _Tag for the evaluation of the derivative on the second dimension._ 
* **struct** [**PolarSplineEvaluator::eval\_type**](structPolarSplineEvaluator_1_1eval__type.md) _Tag for the evaluation of the function._ 
* **struct** [**PolarSplineFEMPoissonLikeSolver::EvalDeriv1DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv1DType.md) _Object storing a value and a value of the derivative of a 1D function._ 
* **struct** [**PolarSplineFEMPoissonLikeSolver::EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) _Object storing a value and a value of the derivatives in each direction of a 2D function._ 
* **struct** [**PolarSplineFEMPoissonLikeSolver::RBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1RBasisSubset.md) 
* **struct** [**PolarSplineFEMPoissonLikeSolver::RCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1RCellDim.md) 
* **struct** [**PolarSplineFEMPoissonLikeSolver::ThetaBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaBasisSubset.md) 
* **struct** [**PolarSplineFEMPoissonLikeSolver::ThetaCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaCellDim.md) 
* **struct** [**PolarSplineMem**](structPolarSplineMem.md) _A structure containing the two FieldMems necessary to define a spline on a set of polar basis splines._ 
* **struct** [**R**](structR.md) _Define non periodic real contravariant_ [_**R**_](structR.md) _dimension._
* **struct** [**R\_cov**](structR__cov.md) _Define non periodic real covariant_ [_**R**_](structR.md) _dimension._
* **struct** [**RefinedDiscreteToCartesianBuilder::Build\_BuilderType**](structRefinedDiscreteToCartesianBuilder_1_1Build__BuilderType.md) 
* **struct** [**RefinedDiscreteToCartesianBuilder::Build\_BuilderType&lt; ddc::SplineBuilder2D&lt; ExecSpace, MemorySpace, BSplinesROriginal, BSplinesThetaOriginal, GridROriginal, GridThetaOriginal, BcLower1, BcUpper1, BcLower2, BcUpper2, Solver, GridROriginal, GridThetaOriginal &gt; &gt;**](structRefinedDiscreteToCartesianBuilder_1_1Build__BuilderType_3_01ddc_1_1SplineBuilder2D_3_01Exe5e607dcd09b0e5dcea3e1a165ab7545e.md) 
* **struct** [**T**](structT.md) _A class which describes the real space in the temporal direction._ 
* **struct** [**Theta**](structTheta.md) _Define periodic real contravariant_ [_**Theta**_](structTheta.md) _dimension._
* **struct** [**Theta\_cov**](structTheta__cov.md) _Define periodic real covariant_ [_**Theta**_](structTheta.md) _dimension._
* **struct** [**Vpar**](structVpar.md) _Define non periodic parallel velocity_  _._
* **struct** [**Vr**](structVr.md) _Define non periodic real_ [_**R**_](structR.md) _velocity dimension._
* **struct** [**Vtheta**](structVtheta.md) _Define periodic real_ [_**Theta**_](structTheta.md) _velocity dimension._
* **struct** [**Vx**](structVx.md) _Define non periodic real_ [_**X**_](structX.md) _velocity dimension._
* **struct** [**Vy**](structVy.md) _Define non periodic real_ [_**Y**_](structY.md) _velocity dimension._
* **struct** [**X**](structX.md) _Define non periodic real_ [_**X**_](structX.md) _dimension._
* **struct** [**X\_pC**](structX__pC.md) _Tag the first non periodic dimension in the pseudo physical domain (pseudo-Cartesian coordinates)._ 
* **struct** [**Y**](structY.md) _Define non periodic real_ [_**Y**_](structY.md) _dimension._
* **struct** [**Y\_pC**](structY__pC.md) _Tag the second non periodic dimension in the pseudo physical domain (pseudo-Cartesian coordinates)._ 
* **struct** [**connectivity\_details::AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md) _A class which helps insert an element into a type sequence._ 
* **struct** [**connectivity\_details::AddToTypeSeq&lt; ToInsert, TypeSeq, BackInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01BackInsert_01_4.md) _Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the back of the type sequence._
* **struct** [**connectivity\_details::AddToTypeSeq&lt; ToInsert, TypeSeq, FrontInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01FrontInsert_01_4.md) _Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the front of the type sequence._
* **struct** [**connectivity\_details::CollectAllGridsOnDim**](structconnectivity__details_1_1CollectAllGridsOnDim.md) _A class which collects all grids along a given dimension in both directions._ 
* **struct** [**connectivity\_details::CollectAllInterfacesOnDim**](structconnectivity__details_1_1CollectAllInterfacesOnDim.md) _A class which collects all grids along a given dimension in both directions._ 
* **struct** [**connectivity\_details::CollectGridsAlongDim**](structconnectivity__details_1_1CollectGridsAlongDim.md) _A class which collects grids along a given dimension on a specified direction from a starting edge._ 
* **struct** [**connectivity\_details::CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md) _Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._
* **struct** [**connectivity\_details::CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins8ee738c554d8fbbf6bab92ba87dd3b80.md) _Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when the grid has already been identified (due to periodicity)._
* **struct** [**connectivity\_details::CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins70deef724c6e45ed62db534c3d9697ec.md) _Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._
* **struct** [**connectivity\_details::CollectInterfacesAlongDim**](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _A class which collects interfaces along a given dimension on a specified direction from a starting edge._ 
* **struct** [**connectivity\_details::CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md) _Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._
* **struct** [**connectivity\_details::CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_36879d7a164b5ac728612e4a981c6d65.md) _Specialisation of_ [_**CollectInterfacesAlongDim**_](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _to stop when the interface has already been identified (due to periodicity)._
* **struct** [**connectivity\_details::CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_ebd86d7b2345baf351562d16964c47d9.md) _Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._
* **struct** [**connectivity\_details::EnforceFirstInterfaceEdge**](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _A class to flip the edges in an interface to ensure that the correct edge comes first._ 
* **struct** [**connectivity\_details::EnforceFirstInterfaceEdge&lt; Interface&lt; Edge2, FirstEdge, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_01Interface_3_01Edge2_00_01FirstEdge_221a02b03250a49af1745b2263467420.md) _Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which needs rearranging._
* **struct** [**connectivity\_details::EnforceFirstInterfaceEdge&lt; Interface&lt; FirstEdge, Edge2, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_01Interface_3_01FirstEdge_00_01Edge2_788676fcb3310ca4c1ec984ff0b4531b.md) _Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which is already correctly arranged._
* **struct** [**connectivity\_details::ExtractPatches**](structconnectivity__details_1_1ExtractPatches.md) _A class to find all the patches used by the various edges._ 
* **struct** [**connectivity\_details::ExtractPatches&lt; ddc::detail::TypeSeq&lt; EdgeType1, EdgeTypes... &gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_01EdgeTypes_8_8_8_01_4_01_4.md) _Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _to iterate recursively over the edge type sequence._
* **struct** [**connectivity\_details::ExtractPatches&lt; ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) _Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _for an empty patch list._
* **struct** [**connectivity\_details::FindInterface**](structconnectivity__details_1_1FindInterface.md) _A class to locate an interface which contains the specified edge._ 
* **struct** [**connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfac6d31b188ee73012ad6c98be99219379f.md) _Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _to iterate recursively over the interface type sequence._
* **struct** [**connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; Edge, OEdge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacd1aa547d7cc4bf022e85928246ab2d07.md) _Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the first interface matches_[_**Edge**_](structEdge.md) _._
* **struct** [**connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacee698732bdf35f06db097afe1714904c.md) _Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the second interface matches_[_**Edge**_](structEdge.md) _._
* **struct** [**connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) _Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for an empty interface list._
* **struct** [**connectivity\_details::FindPatch**](structconnectivity__details_1_1FindPatch.md) _A class to locate a patch which contains the specified grid._ 
* **struct** [**connectivity\_details::FindPatch&lt; Grid1D, ddc::detail::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patch1_00_33770856242f7c5cee1ce419b2efaf64.md) _Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _to iterate recursively over the patch type sequence._
* **struct** [**connectivity\_details::FindPatch&lt; Grid1D, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) _Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _for an empty patch list._
* **struct** [**connectivity\_details::FindPatch&lt; QueryGrid1D, ddc::detail::TypeSeq&lt; Patch&lt; OGrid, QueryGrid1D, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01QueryGrid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patch5f5acd76cfd59a22ebf513823679a320.md) 
* **struct** [**connectivity\_details::FindPatch&lt; QueryGrid1D, ddc::detail::TypeSeq&lt; Patch&lt; QueryGrid1D, OGrid, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01QueryGrid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patchd8fc8921dec760f8fe4c90c2a6947228.md) 
* **struct** [**connectivity\_details::FindRelevantIdxRangeType**](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _A class to find any index range types which contain an index range defined on the provided grid. E.g. Grid1, std::tuple&lt;IdxRange&lt;Grid1, Grid2&gt;, IdxRange&lt;Grid3,Grid4&gt;&gt; will find: ddc::detail::TypeSeq&lt;IdxRange&lt;Grid1, Grid2&gt;&gt;_ 
* **struct** [**connectivity\_details::FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt; IdxRangeHead, IdxRangeTypes... &gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_01Id3b131c802b30082f4412eb4689d6d53b.md) _Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _to iterate recursively over the possible index range types._
* **struct** [**connectivity\_details::FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt;&gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_4_01_4.md) _Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _for an empty list of index range types._
* **struct** [**connectivity\_details::PatchConnection**](structconnectivity__details_1_1PatchConnection.md) _A class which finds all interfaces connected to a given patch._ 
* **struct** [**connectivity\_details::PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01InterfaceType_01_4_01_4.md) _Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an interface list with one element._
* **struct** [**connectivity\_details::PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01Interd9a0a5e7aafe0b71fe7c76720b7c5da6.md) _Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _to iterate recursively over the interface type sequence._
* **struct** [**connectivity\_details::PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) _Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an empty interface list._
* **struct** [**connectivity\_details::SelectRelevantIdxRangeType**](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) _A class to create a type sequence which contains the index range if it can be used to index the grid._ 
* **struct** [**connectivity\_details::SelectRelevantIdxRangeType&lt; QueryGrid1D, IdxRange&lt; IdxRangeGrids... &gt; &gt;**](structconnectivity__details_1_1SelectRelevantIdxRangeType_3_01QueryGrid1D_00_01IdxRange_3_01IdxRangeGrids_8_8_8_01_4_01_4.md) _Specialisation of_ [_**SelectRelevantIdxRangeType**_](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) _to get access to the grids in the index range._
* **struct** [**connectivity\_details::StripOutsideEdges**](structconnectivity__details_1_1StripOutsideEdges.md) _A class which finds all edges which are not_ [_**OutsideEdge**_](structOutsideEdge.md) _types._
* **struct** [**connectivity\_details::StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType_01_4_01_4.md) _Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _for the case with one edge in the list._
* **struct** [**connectivity\_details::StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType1, RemainingEdgeTypes... &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_036e9ce7e4506982efa52c09ca049ae90.md) _Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _to iterate recursively over the edge type sequence._
* **struct** [**connectivity\_details::SwapExtremity**](structconnectivity__details_1_1SwapExtremity.md) _A class to get the opposite edge of a grid line from one of the edges._ 
* **struct** [**connectivity\_details::SwapExtremity&lt; Edge&lt; Patch, Grid1D, BACK &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01BACK_01_4_01_4.md) _Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the back end of a grid line._
* **struct** [**connectivity\_details::SwapExtremity&lt; Edge&lt; Patch, Grid1D, FRONT &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01FRONT_01_4_01_4.md) _Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the front end of a grid line._
* **struct** [**connectivity\_details::ToTuple**](structconnectivity__details_1_1ToTuple.md) _A class to convert a type sequence to a tuple type._ 
* **struct** [**connectivity\_details::ToTuple&lt; ddc::detail::TypeSeq&lt; I... &gt; &gt;**](structconnectivity__details_1_1ToTuple_3_01ddc_1_1detail_1_1TypeSeq_3_01I_8_8_8_01_4_01_4.md) _Specialisation of_ [_**ToTuple**_](structconnectivity__details_1_1ToTuple.md) _for type sequences._
* **struct** [**interpolator\_on\_idx\_range**](structinterpolator__on__idx__range.md) 
* **struct** [**interpolator\_on\_idx\_range&lt; Interp, GridInterp, IdxRange&lt; Grid1D... &gt; &gt;**](structinterpolator__on__idx__range_3_01Interp_00_01GridInterp_00_01IdxRange_3_01Grid1D_8_8_8_01_4_01_4.md) 
* **struct** [**tensor\_tools::GetContravariantDims**](structtensor__tools_1_1GetContravariantDims.md) _A class to get a VectorIndexSet containing only contravariant dimensions._ 
* **struct** [**tensor\_tools::GetContravariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetContravariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) _A class to get a VectorIndexSet containing only contravariant dimensions._ 
* **struct** [**tensor\_tools::GetCovariantDims**](structtensor__tools_1_1GetCovariantDims.md) _A class to get a VectorIndexSet containing only covariant dimensions._ 
* **struct** [**tensor\_tools::GetCovariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetCovariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) _A class to get a VectorIndexSet containing only covariant dimensions._ 
* **struct** [**tensor\_tools::VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) _A class representing a vector index identifier._ 
* **struct** [**tensor\_tools::is\_contravariant\_vector\_index\_set**](structtensor__tools_1_1is__contravariant__vector__index__set.md) 
* **struct** [**tensor\_tools::is\_contravariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__contravariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) _A helper structure to check if all the dimensions in a VectorIndexSet can represent contravariant indices._ 
* **struct** [**tensor\_tools::is\_covariant\_vector\_index\_set**](structtensor__tools_1_1is__covariant__vector__index__set.md) 
* **struct** [**tensor\_tools::is\_covariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__covariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) _A helper structure to check if all the dimensions in a VectorIndexSet can represent covariant indices._ 
* **struct** [**tensor\_tools::vector\_index\_set\_dual**](structtensor__tools_1_1vector__index__set__dual.md) _A helper structure to find a VectorIndexSet describing the covariant indices from a VectorIndexSet describing contravariant indices or to find a VectorIndexSet describing the contravariant indices from a VectorIndexSet describing covariant indices._ 
* **struct** [**tensor\_tools::vector\_index\_set\_dual&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1vector__index__set__dual_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) _The implementation of_ [_**vector\_index\_set\_dual**_](structtensor__tools_1_1vector__index__set__dual.md) _for a VectorIndexSet._
* **class** **IAdvectionV< IdxRangeFdistribu, IdxRangeSpatial >**    
    * **class** [**NullAdvectionVelocity**](classNullAdvectionVelocity.md) _This is a class which imitates a velocity advection. It inherits from IAdvectionV and can be used as an advection operator but does not actually modify the distribution function. This can be useful for debugging purposes._ 
* **class** **std::conditional_t< BsplineOnUniformCellsMu, ddc::UniformBSplines< Mu, BSDegreeMu >, ddc::NonUniformBSplines< Mu, BSDegreeMu > >**    
    * **struct** [**BSplinesMu**](structBSplinesMu.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsR, ddc::UniformBSplines< R, BSDegreeR >, ddc::NonUniformBSplines< R, BSDegreeR > >**    
    * **struct** [**BSplinesR**](structBSplinesR.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsTheta, ddc::UniformBSplines< Theta, BSDegreeTheta >, ddc::NonUniformBSplines< Theta, BSDegreeTheta > >**    
    * **struct** [**BSplinesTheta**](structBSplinesTheta.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsVpar, ddc::UniformBSplines< Vpar, BSDegreeVpar >, ddc::NonUniformBSplines< Vpar, BSDegreeVpar > >**    
    * **struct** [**BSplinesVpar**](structBSplinesVpar.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsVx, ddc::UniformBSplines< Vx, BSDegreeVx >, ddc::NonUniformBSplines< Vx, BSDegreeVx > >**    
    * **struct** [**BSplinesVx**](structBSplinesVx.md) 
    * **struct** [**BSplinesVx**](structBSplinesVx.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsVy, ddc::UniformBSplines< Vy, BSDegreeVy >, ddc::NonUniformBSplines< Vy, BSDegreeVy > >**    
    * **struct** [**BSplinesVy**](structBSplinesVy.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsX, ddc::UniformBSplines< X, BSDegreeX >, ddc::NonUniformBSplines< X, BSDegreeX > >**    
    * **struct** [**BSplinesX**](structBSplinesX.md) 
    * **struct** [**BSplinesX**](structBSplinesX.md) 
    * **struct** [**BSplinesX**](structBSplinesX.md) 
* **class** **std::conditional_t< BsplineOnUniformCellsY, ddc::UniformBSplines< Y, BSDegreeY >, ddc::NonUniformBSplines< Y, BSDegreeY > >**    
    * **struct** [**BSplinesY**](structBSplinesY.md) 
    * **struct** [**BSplinesY**](structBSplinesY.md) 
* **class** **std::conditional_t< uniform_edge_v, UniformGridBase< Vx >, NonUniformGridBase< Vx > >**    
    * **struct** [**CollisionsIntra::GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) 
    * **struct** [**CollisionsIntra::GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) 
* **class** **NonUniformGridBase< PDEDim >**    
    * **struct** [**FEM1DPoissonSolver::GridPDEDimQ**](structFEM1DPoissonSolver_1_1GridPDEDimQ.md) _The grid of quadrature points along the PDEDim direction._ 
* **class** **ddc::NonUniformBSplines< PDEDim, InputBSplines::degree()>**    
    * **struct** [**FEM1DPoissonSolver::HiddenFEMBSplines**](structFEM1DPoissonSolver_1_1HiddenFEMBSplines.md) 
* **class** **ddc::PeriodicSampling< ddc::Fourier< Dim > >**    
    * **struct** [**FFTPoissonSolver&lt; IdxRange&lt; GridPDEDim1D... &gt;, IdxRangeFull, ExecSpace, LayoutSpace &gt;::GridFourier**](structFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace2aeecfe91d464f5738599cc105fb6087.md) 
* **class** **SplineInterpPointsMu::interpolation_discrete_dimension_type**    
    * **struct** [**GridMu**](structGridMu.md) 
* **class** **SplineInterpPointsR::interpolation_discrete_dimension_type**    
    * **struct** [**GridR**](structGridR.md) 
* **class** **SplineInterpPointsTheta::interpolation_discrete_dimension_type**    
    * **struct** [**GridTheta**](structGridTheta.md) 
* **class** **SplineInterpPointsVpar::interpolation_discrete_dimension_type**    
    * **struct** [**GridVpar**](structGridVpar.md) 
* **class** **SplineInterpPointsVx::interpolation_discrete_dimension_type**    
    * **struct** [**GridVx**](structGridVx.md) 
    * **struct** [**GridVx**](structGridVx.md) 
* **class** **SplineInterpPointsVy::interpolation_discrete_dimension_type**    
    * **struct** [**GridVy**](structGridVy.md) 
* **class** **SplineInterpPointsX::interpolation_discrete_dimension_type**    
    * **struct** [**GridX**](structGridX.md) 
    * **struct** [**GridX**](structGridX.md) 
    * **struct** [**GridX**](structGridX.md) 
* **class** **SplineInterpPointsY::interpolation_discrete_dimension_type**    
    * **struct** [**GridY**](structGridY.md) 
    * **struct** [**GridY**](structGridY.md) 
* **class** **NonUniformGridBase< R >**    
    * **struct** [**PolarSplineFEMPoissonLikeSolver::QDimRMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimRMesh.md) _Tag the first dimension for the quadrature mesh._ 
* **class** **NonUniformGridBase< Theta >**    
    * **struct** [**PolarSplineFEMPoissonLikeSolver::QDimThetaMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimThetaMesh.md) _Tag the second dimension for the quadrature mesh._ 
* **class** **std::conditional_t< BSplinesROriginal::is_uniform(), ddc::UniformBSplines< R, BSplinesROriginal::degree()>, ddc::NonUniformBSplines< R, BSplinesROriginal::degree()> >**    
    * **struct** [**RefinedDiscreteToCartesianBuilder::BSplinesRRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesRRefined.md) _The type of the radial B-splines on which the new mapping will be defined._ 
* **class** **std::conditional_t< BSplinesThetaOriginal::is_uniform(), ddc::UniformBSplines< Theta, BSplinesThetaOriginal::degree()>, ddc::NonUniformBSplines< Theta, BSplinesThetaOriginal::degree()> >**    
    * **struct** [**RefinedDiscreteToCartesianBuilder::BSplinesThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesThetaRefined.md) _The type of the poloidal B-splines on which the new mapping will be defined._ 
* **class** **GrevillePointsR::interpolation_discrete_dimension_type**    
    * **struct** [**RefinedDiscreteToCartesianBuilder::GridRRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridRRefined.md) _The type of the grid of radial points on which the new mapping will be defined._ 
* **class** **GrevillePointsTheta::interpolation_discrete_dimension_type**    
    * **struct** [**RefinedDiscreteToCartesianBuilder::GridThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridThetaRefined.md) _The type of the grid of poloidal points on which the new mapping will be defined._ 
* **class** **std::false_type**    
    * **struct** [**ddcHelper::is\_non\_uniform\_interpolation\_points**](structddcHelper_1_1is__non__uniform__interpolation__points.md) 
    * **struct** [**is\_onion\_patch\_locator**](structis__onion__patch__locator.md) _Struct to identify if the patch locator is adapted to onion geometry._ 
    * **struct** [**is\_subidx\_range\_collection**](structis__subidx__range__collection.md) 
    * **struct** [**tensor\_tools::is\_vector\_index\_set**](structtensor__tools_1_1is__vector__index__set.md) _A helper structure to recognise a VectorIndexSet type._ 
* **class** **std::true_type**    
    * **struct** [**ddcHelper::is\_non\_uniform\_interpolation\_points&lt; NonUniformInterpolationPoints&lt; BSplines, BcXmin, BcXmax &gt; &gt;**](structddcHelper_1_1is__non__uniform__interpolation__points_3_01NonUniformInterpolationPoints_3_047d1c8570873e3c052e2e394afcf9270.md) 
    * **struct** [**is\_onion\_patch\_locator&lt; OnionPatchLocator&lt; MultipatchIdxRanges, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace &gt; &gt;**](structis__onion__patch__locator_3_01OnionPatchLocator_3_01MultipatchIdxRanges_00_01LogicalToPhys15c96379834346672a2b2d644897e91f.md) 
    * **struct** [**is\_subidx\_range\_collection&lt; IdxRangeSlice&lt; Tags... &gt; &gt;**](structis__subidx__range__collection_3_01IdxRangeSlice_3_01Tags_8_8_8_01_4_01_4.md) 
    * **struct** [**tensor\_tools::is\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) 

