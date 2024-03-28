#pragma once

#include <ddc/ddc.hpp>

#include <KOLIOP/koliop.h>

#include "assert.hpp"
#include "geometry.hpp"

/**
 * @brief Header Guard for collision operator
 */
class CollisionsGuard
{
public:
    /**
      * @brief
      * @param[in] a_host_thread_count When compiled for CPU, this represents
      * the number of OpenMP thread the operator should use.
      * @param[in] a_device_index The device the operator should use, you almost
      * always want to use the first and only one visible.
     */
    CollisionsGuard(std::size_t a_host_thread_count, std::size_t a_device_index = 0)
    {
        if (::koliop_Initialize(false, a_host_thread_count, a_device_index, true)
            != KOLIOP_STATUS_SUCCESS) {
            GSLX_ASSERT(false);
        }
    }

    ~CollisionsGuard()
    {
        if (::koliop_Deinitialize(false) != KOLIOP_STATUS_SUCCESS) {
            GSLX_ASSERT(false);
        }
    }
};


/**
 * @brief A class which computes the collision operator in (vpar,mu)
 */
class Collisions /* : public IRightHandSide */
{
public:
    /**
     * @brief Create instance of Collisions class
     *
     * @param[in] dom_sp_tor3D_v2D
     *      domain (species, 3D space, 2D velocity space) on which the collision operator acts
     *      ATTENTION it must contain the all domain in species and 2D velocity.
     * @param[in] coeff_intdmu
     *      quadrature coefficient in mu direction
     * @param[in] coeff_intdvpar
     *      quadrature coefficient in vpar direction
     * @param[in] nustar0_r
     *      radial profile of collisionality
     */
    Collisions(
            DDomSpTor3DV2D dom_sp_tor3D_v2D,
            device_t<DViewMu> coeff_intdmu,
            device_t<DViewVpar> coeff_intdvpar,
            device_t<DViewTor1> nustar0_r);
    ~Collisions();

    /**
     * @brief Apply the collision operator to the distribution functions of all species on all species
     *
     * @param[inout] all_f_distribution
     *      All distribution functions
     * @param[in] deltat_coll
     *      Collision time step
     */
    void operator()(device_t<DSpanSpTor3DV2D> all_f_distribution, double deltat_coll) const
            /* override */;

protected:
    /**
     * Opaque type representing the operator (due to the C interface)
    */
    ::koliop_Operator m_operator_handle;

    /**
     * MDL is a helper type which avoids repetition of a lengthy type name.
     * In the future this type may be removed once the objects which use these
     * types are defined elsewhere with DDC objects.
     * MDL stands for:
     * M : Managed - This indicates that the type owns its data and is responsible for deallocating it.
     * D : Device - The memory is saved on GPU (when compiled for GPU).
     * L : Left memory layout - The multi-indices are strided left-to-right (meaning indexing matches
     * 			the standard indexing in Fortran not C)
     */
    template <typename Extent>
    using MDL = Kokkos::View<
            Extent,
            Kokkos::LayoutLeft,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::MemoryTraits<Kokkos::Restrict>>;

    // NOTE: Some of these arrays should come from a parent class or manager.
    // They resides in this class while we wait for their implementation.
    /// Combinatory (6x6) matrix computed only one times at initialisation. Rk: 6 = 2*(Npolmax-1) + 1 + 1
    MDL<double[6][6]> m_comb_mat;
    /// Normalized masses for all species
    MDL<double*> m_hat_As; // DFieldSp
    /// Normalized charges for all species
    MDL<double*> m_hat_Zs;
    /// Mesh points in the radial direction
    // TODO: See if we need m_rg ?
    MDL<double*> m_rg;
    /// Radial safety factor profile
    MDL<double*> m_q_r;
    /// Mask used to avoid to apply collision in certain region
    // TODO: Understand this mask
    MDL<double*> m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    MDL<double**> m_mask_LIM;
    /// B norm in (r,theta)
    // TODO: Attention this must be 3D --> transfer it in a 1D array ?
    MDL<double**> m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // TODO : Must be 5D for full 3D geometry
    MDL<double****> m_Bstar_s;
    /// grid in mu direction
    MDL<double*> m_mug;
    /// grid in vpar direction
    MDL<double*> m_vparg;
};
