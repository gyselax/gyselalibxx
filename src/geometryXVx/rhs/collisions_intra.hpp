// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/matrix_banded.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

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
    static constexpr bool uniform_edge_v = ddc::is_uniform_point_sampling_v<GridVx>;

public:
    /**
     * A conditional type representing either a uniform or a non-uniform ghosted vx mesh. 
     */
    struct GhostedVx
        : std::conditional_t<uniform_edge_v, UniformGridBase<Vx>, NonUniformGridBase<Vx>>
    {
    };

    /**
     * A conditional type representing either a uniform or a non-uniform ghosted staggered vx mesh. 
     */
    struct GhostedVxStaggered
        : std::conditional_t<uniform_edge_v, UniformGridBase<Vx>, NonUniformGridBase<Vx>>
    {
    };

    /**
     * A type representing a mesh for species, space and ghosted vx mesh. 
     */
    using IdxRangeSpXVx_ghosted = IdxRange<Species, GridX, GhostedVx>;

    /**
     * A type representing a mesh for species, space and ghosted staggered vx mesh. 
     */
    using IdxRangeSpXVx_ghosted_staggered = IdxRange<Species, GridX, GhostedVxStaggered>;

    /**
     * A type representing a ghosted vx index. 
     */
    using IndexVx_ghosted = Idx<GhostedVx>;

    /**
     * A type representing a ghosted staggered vx index. 
     */
    using IndexVx_ghosted_staggered = Idx<GhostedVxStaggered>;

    /**
     * A type representing a species, space and ghosted vx index. 
     */
    using IndexSpXVx_ghosted = Idx<Species, GridX, GhostedVx>;

    /**
     * A type representing a species, space and ghosted staggered vx index. 
     */
    using IndexSpXVx_ghosted_staggered = Idx<Species, GridX, GhostedVxStaggered>;


private:
    template <class TargetDim>
    KOKKOS_FUNCTION static Idx<TargetDim> to_index(Idx<GridVx> const& index);

    template <class VDim>
    std::enable_if_t<!ddc::is_uniform_point_sampling_v<VDim>>
    build_ghosted_staggered_vx_point_sampling(IdxRange<VDim> const& dom);

    template <class VDim>
    std::enable_if_t<ddc::is_uniform_point_sampling_v<VDim>>
    build_ghosted_staggered_vx_point_sampling(IdxRange<VDim> const& dom);

    double m_nustar0;
    double m_fthresh;
    DFieldMemSpX m_nustar_profile_alloc;
    DFieldSpX m_nustar_profile;

    IdxRange<GhostedVx> m_gridvx_ghosted;
    IdxRange<GhostedVxStaggered> m_gridvx_ghosted_staggered;

    IdxRangeSpXVx_ghosted m_mesh_ghosted;
    IdxRangeSpXVx_ghosted_staggered m_mesh_ghosted_staggered;

public:
    /**
     * @brief The constructor for the operator.
     *
     * @param[in] mesh The index range on which the operator will act.
     * @param[in] nustar0 The normalized collisionality.
     */
    CollisionsIntra(IdxRangeSpXVx const& mesh, double nustar0);

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
     * @return A field referencing the distribution function passed as argument.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;

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
    IdxRange<GhostedVx> const& get_gridvx_ghosted() const;

    /**
     * @brief Get the ghosted and staggered vx mesh used for computing finite differences centered derivatives.
     *
     * @return The ghosted and staggered vx mesh.
     */
    IdxRange<GhostedVxStaggered> const& get_gridvx_ghosted_staggered() const;

    /**
     * @brief Get a mesh containing the species, spatial and the ghosted vx mesh.
     *
     * @return The species, spatial, and ghosted vx mesh.
     */
    IdxRange<Species, GridX, GhostedVx> const& get_mesh_ghosted() const;

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
            DFieldSpXVx RR,
            DConstFieldSpXVx AA,
            DConstFieldSpXVx BB,
            DConstFieldSpXVx CC,
            DConstFieldSpXVx allfdistribu,
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
            DFieldSpXVx AA,
            DFieldSpXVx BB,
            DFieldSpXVx CC,
            Field<double, IdxRangeSpXVx_ghosted> Dcoll,
            Field<double, IdxRangeSpXVx_ghosted_staggered> Dcoll_staggered,
            Field<double, IdxRangeSpXVx_ghosted> Nucoll,
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
            host_t<DConstFieldVx> AA,
            host_t<DConstFieldVx> BB,
            host_t<DConstFieldVx> CC) const;
};
