#include <KOLIOP/koliop.h>

#include "assert.hpp"
#include "collision_operator.hpp"

// NOTE: Offset of in the r dimension. Used in the Fortran code.
static inline constexpr std::size_t local_idx_range_r_offset = 0;
// NOTE: Offset of in the theta dimension. Used in the Fortran code.
static inline constexpr std::size_t local_idx_range_theta_offset = 0;

// NOTE: For now, these are constant. When more physics get implemented, pass
// them to the Collisions operator's constructor.
// NOTE: Enable masking.
static inline constexpr std::int8_t is_sol_enabled = false;
static inline constexpr std::int8_t is_limiter_mask_enabled = true;

namespace detail {

::koliop_Operator do_operator_initialisation(
        std::size_t mu_extent,
        std::size_t vpar_extent,
        std::size_t tor1_extent,
        std::size_t tor2_extent,
        std::size_t tor3_extent,
        std::size_t species_extent,
        std::int8_t collision_interspecies,
        std::int64_t ir_sol_separatrix,
        double const* hat_As,
        double const* hat_Zs,
        double const* mug,
        double const* vparg,
        double const* coeff_intdmu,
        double const* coeff_intdvpar,
        double const* coeff_AD_r,
        double const* mask_buffer_r,
        double const* mask_limiter,
        double const* B_norm,
        double const* Bstar_s)
{
    ::koliop_Operator operator_handle;

    std::size_t const r_extent = tor2_extent;
    std::size_t const theta_extent = tor1_extent;

    if (::koliop_Create(
                &operator_handle,
                mu_extent,
                vpar_extent,
                tor1_extent,
                tor2_extent,
                tor3_extent,
                species_extent,
                r_extent, // NOTE: No MPI index range decomposition.
                theta_extent, // NOTE: No MPI index range decomposition.
                local_idx_range_r_offset,
                local_idx_range_theta_offset,
                collision_interspecies,
                is_sol_enabled,
                is_limiter_mask_enabled,
                ir_sol_separatrix,
                hat_As,
                hat_Zs,
                mug,
                vparg,
                coeff_intdmu,
                coeff_intdvpar,
                coeff_AD_r,
                mask_buffer_r,
                mask_limiter,
                B_norm,
                Bstar_s)
        != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }

    return operator_handle;
}

void do_operator_deinitialisation(::koliop_Operator operator_handle)
{
    if (::koliop_Release(operator_handle) != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }
}

} // namespace detail
