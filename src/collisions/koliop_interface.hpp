#pragma once

#include <ddc/ddc.hpp>

#include <KOLIOP/koliop.h>

/**
 * @brief A namespace containing functions which call the koliop library
 */
namespace koliop_interface {
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

/**
 * @brief Initialise the collision operator defined in koliop.
 *
 * @param[in] the_mu_extent The number of points in the magnetic moment dimension.
 * @param[in] the_vpar_extent The number of points in the parallel velocity dimension.
 * @param[in] the_r_extent The number of points in the radial dimension.
 * @param[in] the_theta_extent The number of points in the poloidal dimension.
 * @param[in] the_phi_extent The number of points in the toroidal dimension.
 * @param[in] the_species_extent The number of species.
 * @param[in] ir_SOL_separatrix The index in the radial dimension where the SOL separatrix is found.
 * @param[in] mug The coordinates of the grid of magnetic moments.
 * @param[in] vparg The coordinates of the grid of parallel velocities.
 * @param[in] coeff_intdmu The quadrature coefficients to integrate on the grid of magnetic moments.
 * @param[in] coeff_intdvpar The quadrature coefficients to integrate on the grid of parallel velocities.
 * @param[in] nustar The radial profile of the collisionality.
 * @param[in] comb_mat The combinatory (6x6) matrix.
 * @param[in] hat_As The normalized masses for all species.
 * @param[in] hat_Zs The normalized charges for all species.
 * @param[in] rg The coordinates of the grid of radial points.
 * @param[in] safety_factor The radial profile of the safety factor.
 * @param[in] mask_buffer_r The mask used to avoid applying collisions in certain radial region.
 * @param[in] mask_LIM The limiter mask defined over the polar slice.
 * @param[in] B_norm The norm of the magnetic field defined over the polar slice.
 * @param[in] Bstar_s Bstar defined on the coordinates (species,r,theta,vpar).
 *
 * @returns The handle which identifies the operator.
 */
::koliop_Operator DoOperatorInitialization(
        std::size_t the_mu_extent,
        std::size_t the_vpar_extent,
        std::size_t the_r_extent,
        std::size_t the_theta_extent,
        std::size_t the_phi_extent,
        std::size_t the_species_extent,
        std::int64_t ir_SOL_separatrix,
        const double* mug,
        const double* vparg,
        const double* coeff_intdmu,
        const double* coeff_intdvpar,
        const double* nustar,
        const double* comb_mat,
        const double* hat_As,
        const double* hat_Zs,
        const double* rg,
        const double* safety_factor,
        const double* mask_buffer_r,
        const double* mask_LIM,
        const double* B_norm,
        const double* Bstar_s);

/**
 * @brief Calculate the combinatory (6x6) matrix.
 *
 * @param[out] comb_mat The combinatory (6x6) matrix.
 */
void DoCombMatComputation(MDL<double[6][6]>& comb_mat);


/**
 * @brief Destructor for koliop.
 *
 * @param the_operator_handle The handle which identifies the operator.
 */
void DoOperatorDeinitialization(::koliop_Operator the_operator_handle);
}; // namespace koliop_interface
