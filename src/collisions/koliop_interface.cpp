#include <atomic>
#include <cassert>
#include <iostream>

#include <KOLIOP/koliop.h>

#include "assert.hpp"
#include "koliop_interface.hpp"

// NOTE: For now, these are constant. When more physics get implemented, pass
// them via the Collisions' constructor.
static inline constexpr std::size_t Npolmax = 3;
static inline constexpr std::size_t index_max = 2 * (Npolmax - 1);
// pSpecies_id represents the species processed by a rank.
//  It is fixed to 0 but is no longer used since the 2024/03/01 modification in the Fortran version.
//  After this date and on both the C++ and Fortran versions, the value is meaningless as a rank contains all
//  the species for a given global simulation space sub volume.
// [TODO] Delete this input (both for Fortran and C++ versions)
static inline constexpr std::int64_t pSpecies_id = 0;
static inline constexpr std::int64_t mu_id
        = 0; // [TODO] Same than for pSpecies_id: must be deleted from the list of input
static inline constexpr std::size_t the_local_idx_range_r_offset = 0;
static inline constexpr std::size_t the_local_idx_range_theta_offset = 0;
static inline constexpr std::int8_t SOL = false;
static inline constexpr std::int8_t LIM = false;

namespace koliop_interface {

::koliop_Operator DoOperatorInitialization(
        std::size_t the_mu_extent,
        std::size_t the_vpar_extent,
        std::size_t the_r_extent,
        std::size_t the_theta_extent,
        std::size_t the_phi_extent,
        std::size_t the_species_extent,
        std::int8_t collision_interspecies,
        std::int64_t ir_SOL_separatrix,
        const double* hat_As,
        const double* hat_Zs,
        const double* mug,
        const double* vparg,
        const double* coeff_intdmu,
        const double* coeff_intdvpar,
        const double* coeff_AD_r,
        const double* mask_buffer_r,
        const double* mask_LIM,
        const double* B_norm,
        const double* Bstar_s)
{
    ::koliop_Operator the_operator_handle;

    if (::koliop_Create(
                &the_operator_handle,
                the_mu_extent,
                the_vpar_extent,
                the_r_extent,
                the_theta_extent,
                the_phi_extent,
                the_species_extent,
                the_r_extent, // NOTE: No MPI index range decomposition.
                the_theta_extent, // NOTE: No MPI index range decomposition.
                the_local_idx_range_r_offset,
                the_local_idx_range_theta_offset,
                collision_interspecies,
                SOL,
                LIM,
                ir_SOL_separatrix,
                hat_As,
                hat_Zs,
                mug,
                vparg,
                coeff_intdmu,
                coeff_intdvpar,
                coeff_AD_r,
                mask_buffer_r,
                mask_LIM,
                B_norm,
                Bstar_s)
        != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }

    return the_operator_handle;
}

void DoOperatorDeinitialization(::koliop_Operator the_operator_handle)
{
    if (::koliop_Release(the_operator_handle) != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }
}

/**
 * UHL is a helper type which avoids repetition of a lengthy type name.
 * In the future this type may be removed once the objects which use these
 * types are defined elsewhere with DDC objects.
 * UHL stands for:
 * U : Unmanaged - This indicates that the type does not own its data and will not deallocate it.
 * H : Host - The memory is saved on CPU.
 * L : Left memory layout - The multi-indices are strided left-to-right (meaning indexing matches
 * 			the standard indexing in Fortran not C)
 */
template <typename Extent>
using UHL = Kokkos::View<
        Extent,
        Kokkos::LayoutLeft,
        Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>>;

void DoCombMatComputation(MDL<double[6][6]>& comb_mat)
{
    // NOTE: In layout left.
    // Precompute: n! /(k! * (n-k)!)

    // clang-format off
    static constexpr double kCombMat[6 * 6]
            = {1.0, 0.0, 0.0,   0.0, 0.0, 0.0,
               1.0, 1.0, 0.0,   0.0, 0.0, 0.0,
               1.0, 2.0, 1.0,   0.0, 0.0, 0.0,
               1.0, 3.0, 3.0,   1.0, 0.0, 0.0,
               1.0, 4.0, 6.0,   4.0, 1.0, 0.0,
               1.0, 5.0, 10.0, 10.0, 5.0, 1.0};
    // clang-format on

    Kokkos::deep_copy(comb_mat, UHL<const double**>(kCombMat, 6, 6));
}

} // namespace koliop_interface
