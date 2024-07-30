#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/matrix_banded.hpp>

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

/**
 * @brief Class describing the intra-species collision operator
 *  
 * The intra-species collision operator can be written as @f$ C_{aa} = d/dv (Dcoll df_a/dv - Nucoll f_a) @f$
 * Where Dcoll and Nucoll are the collision operator diffusion and advection coefficients 
 * that depend on the absolute value of the velocity and possibly on space. f_a is the distribution 
 * function. The Nucoll and Dcoll coefficients are adjusted at each time and spatial position 
 * so that the collision operator conserves the density, momentum and energy of the 
 * considered species it acts on.
 * 
 * The evolution equation for the collisions @f$ df_a/dt = C_{aa} @f$ is solved using a Crank-Nicolson 
 * finite difference scheme adapted for non-uniform meshes. The derivatives are computed using
 * centered second-order finite differences. To use the same derivatives formula at the edges of 
 * the vx mesh, we introduce a ghosted vx mesh with additional points at the edges. To further 
 * improve the accuracy of the computations of the derivatives, we also introduce another vx 
 * mesh that is staggered with respect to the initial mesh. The points of the staggered mesh 
 * lie at the middle of the cells of the initial vx mesh (where a cell is defined by two 
 * adjacent points). The Crank-Nicolson scheme leads to the formulation of a linear system 
 * that needs to be resolved at each spatial position of the simulation box. Note that this linear 
 * system depends on the considered spatial position. 
 * 
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/collisions_intra_inter.pdf).
 */
class CollisionsIntra : public IRightHandSide
{
private:
    static constexpr bool uniform_edge_v = ddc::is_uniform_point_sampling_v<IDimVx>;

public:
    /**
     * A conditional type representing either a uniform or a non-uniform ghosted vx mesh. 
     */
    struct GhostedVx
        : std::conditional_t<
                  uniform_edge_v,
                  ddc::UniformPointSampling<RDimVx>,
                  ddc::NonUniformPointSampling<RDimVx>>
    {
    };

    /**
     * A conditional type representing either a uniform or a non-uniform ghosted staggered vx mesh. 
     */
    struct GhostedVxStaggered
        : std::conditional_t<
                  uniform_edge_v,
                  ddc::UniformPointSampling<RDimVx>,
                  ddc::NonUniformPointSampling<RDimVx>>
    {
    };

    /**
     * A type representing a mesh for species, space and ghosted vx mesh. 
     */
    using IDomainSpXVx_ghosted = ddc::DiscreteDomain<IDimSp, IDimX, GhostedVx>;

    /**
     * A type representing a mesh for species, space and ghosted staggered vx mesh. 
     */
    using IDomainSpXVx_ghosted_staggered = ddc::DiscreteDomain<IDimSp, IDimX, GhostedVxStaggered>;

    /**
     * A type representing a ghosted vx index. 
     */
    using IndexVx_ghosted = ddc::DiscreteElement<GhostedVx>;

    /**
     * A type representing a ghosted staggered vx index. 
     */
    using IndexVx_ghosted_staggered = ddc::DiscreteElement<GhostedVxStaggered>;

    /**
     * A type representing a species, space and ghosted vx index. 
     */
    using IndexSpXVx_ghosted = ddc::DiscreteElement<IDimSp, IDimX, GhostedVx>;

    /**
     * A type representing a species, space and ghosted staggered vx index. 
     */
    using IndexSpXVx_ghosted_staggered = ddc::DiscreteElement<IDimSp, IDimX, GhostedVxStaggered>;


private:
    template <class TargetDim>
    KOKKOS_FUNCTION static ddc::DiscreteElement<TargetDim> to_index(
            ddc::DiscreteElement<IDimVx> const& index);

    template <class VDim>
    std::enable_if_t<!ddc::is_uniform_point_sampling_v<VDim>>
    build_ghosted_staggered_vx_point_sampling(ddc::DiscreteDomain<VDim> const& dom);

    template <class VDim>
    std::enable_if_t<ddc::is_uniform_point_sampling_v<VDim>>
    build_ghosted_staggered_vx_point_sampling(ddc::DiscreteDomain<VDim> const& dom);

    double m_nustar0;
    double m_fthresh;
    DFieldSpX m_nustar_profile_alloc;
    DSpanSpX m_nustar_profile;

    ddc::DiscreteDomain<GhostedVx> m_gridvx_ghosted;
    ddc::DiscreteDomain<GhostedVxStaggered> m_gridvx_ghosted_staggered;

    IDomainSpXVx_ghosted m_mesh_ghosted;
    IDomainSpXVx_ghosted_staggered m_mesh_ghosted_staggered;

public:
    /**
     * @brief The constructor for the operator.
     *
     * @param[in] mesh The domain on which the operator will act.
     * @param[in] nustar0 The normalized collisionality.
     */
    CollisionsIntra(IDomainSpXVx const& mesh, double nustar0);

    ~CollisionsIntra() = default;

    /**
     * @brief Update the distribution function for intra-species collision.
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following collisions within various species.
     * This operator only handles collisions between particles of the same
     * species.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step over which the collisions occur.
     *
     * @return A span referencing the distribution function passed as argument.
     */
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

    /**
     * @brief Get the collision coefficient.
     *
     * @return The collisionality coefficient.
     */
    double get_nustar0() const;

    /**
     * @brief Get the ghosted vx mesh used for computing finite differences centered derivatives.
     *
     * @return The ghosted vx mesh.
     */
    ddc::DiscreteDomain<GhostedVx> const& get_gridvx_ghosted() const;

    /**
     * @brief Get the ghosted and staggered vx mesh used for computing finite differences centered derivatives.
     *
     * @return The ghosted and staggered vx mesh.
     */
    ddc::DiscreteDomain<GhostedVxStaggered> const& get_gridvx_ghosted_staggered() const;

    /**
     * @brief Get a mesh containing the species, spatial and the ghosted vx mesh.
     *
     * @return The species, spatial, and ghosted vx mesh.
     */
    ddc::DiscreteDomain<IDimSp, IDimX, GhostedVx> const& get_mesh_ghosted() const;

    /**
     * @brief Compute the right-hand-side of the collision operator linear system.
     * @param[inout] RR A vector representing the right-hand-side of the system.
     * @param[in] AA A vector representing the lower diagonal of the matrix of the linear system.
     * @param[in] BB A vector representing the diagonal of the matrix of the linear system.
     * @param[in] CC A vector representing the upper diagonal of the matrix of the linear system.
     * @param[in] allfdistribu The distribution function. 
     * @param[in] fthresh A constant value used for imposing Dirichlet boundary condition to solve the linear system.  
     */
    void compute_rhs_vector(
            DSpanSpXVx RR,
            DViewSpXVx AA,
            DViewSpXVx BB,
            DViewSpXVx CC,
            DViewSpXVx allfdistribu,
            double fthresh) const;
    /**
     * @brief Compute the coefficients of the tridiagonal matrix of the collision operator linear system.
     * @param[inout] AA A vector representing the lower diagonal of the matrix of the linear system.
     * @param[inout] BB A vector representing the diagonal of the matrix of the linear system.
     * @param[inout] CC A vector representing the upper diagonal of the matrix of the linear system.
     * @param[in] Dcoll The velocity-dependent diffusion coefficient of the collision operaror. 
     * @param[in] Dcoll_staggered The velocity-dependent diffusion coefficient of the collision operaror computed on the staggered vx mesh. 
     * @param[in] Nucoll The velocity-dependent advection coefficient of the collision operaror. 
     * @param[in] deltat The time step.
     */
    void compute_matrix_coeff(
            DSpanSpXVx AA,
            DSpanSpXVx BB,
            DSpanSpXVx CC,
            device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Dcoll,
            device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted_staggered>> Dcoll_staggered,
            device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Nucoll,
            double deltat) const;

    /**
     * @brief Fill the tridiagonal matrix of the collision operator linear system with the required coefficients.
     * @param[inout] matrix A banded (tridiagonal) matrix.
     * @param[in] AA A vector containing the lower diagonal coefficients of the matrix.
     * @param[in] BB A vector containing the diagonal coefficients of the matrix.
     * @param[in] CC A vector containing the upper diagonal coefficients of the matrix.
     */
    void fill_matrix_with_coeff(
            Matrix_Banded& matrix,
            host_t<DViewVx> AA,
            host_t<DViewVx> BB,
            host_t<DViewVx> CC) const;
};
