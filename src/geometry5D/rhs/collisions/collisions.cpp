#include <atomic>
#include <cassert>
#include <iostream>

#include <KOLIOP/koliop.h>

#include "assert.hpp"
#include "collisions.hpp"

// NOTE: For now, these are constant. When more physics get implemented, pass
// them via the Collisions' constructor.
static inline constexpr std::size_t Npolmax = 3;
static inline constexpr std::size_t index_max = 2 * (Npolmax - 1);
static inline constexpr std::int64_t pSpecies_id = 0;
static inline constexpr std::int64_t mu_id = 0;
static inline constexpr std::size_t the_local_domain_r_offset = 0;
static inline constexpr std::size_t the_local_domain_theta_offset = 0;
static inline constexpr double R0 = 1.0;
static inline constexpr std::int8_t collision_interspecies = false;
static inline constexpr std::int8_t SOL = false;
static inline constexpr std::int8_t LIM = false;

namespace detail {

static inline ::koliop_Operator DoOperatorInitialization(
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
        const double* nustar0_r,
        const double* comb_mat,
        const double* hat_As,
        const double* hat_Zs,
        const double* rg,
        const double* q_r,
        const double* mask_buffer_r,
        const double* mask_LIM,
        const double* B_norm,
        const double* Bstar_s)
{
    ::koliop_Operator the_operator_handle;

    if (::koliop_Create(
                &the_operator_handle,
                Npolmax,
                index_max,
                pSpecies_id,
                mu_id,
                the_mu_extent,
                the_vpar_extent,
                the_r_extent,
                the_theta_extent,
                the_phi_extent,
                the_species_extent,
                the_r_extent, // NOTE: No MPI domain decomposition.
                the_theta_extent, // NOTE: No MPI domain decomposition.
                the_local_domain_r_offset,
                the_local_domain_theta_offset,
                R0,
                collision_interspecies,
                SOL,
                LIM,
                ir_SOL_separatrix,
                comb_mat,
                hat_As,
                hat_Zs,
                mug,
                vparg,
                coeff_intdmu,
                coeff_intdvpar,
                rg,
                nustar0_r,
                q_r,
                mask_buffer_r,
                mask_LIM,
                B_norm,
                Bstar_s)
        != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }

    return the_operator_handle;
}

static inline void DoOperatorDeinitialization(::koliop_Operator the_operator_handle)
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
        Kokkos::DefaultHostExecutionSpace::memory_space,
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>>;

template <typename View2D>
static inline void DoCombMatComputation(View2D& comb_mat)
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

template <typename Span1D, typename CoordinateSampler>
static inline void DumpCoordinates(
        Span1D& dump_coord,
        ddc::DiscreteDomain<CoordinateSampler> sampler)
{
    ddc::ChunkSpan ddc_dump_coord(dump_coord, sampler);
    ddc::parallel_for_each(
            typename Span1D::execution_space {},
            sampler,
            KOKKOS_LAMBDA(ddc::DiscreteElement<CoordinateSampler> i) {
                ddc_dump_coord(i) = ddc::coordinate(i);
            });
}

} // namespace detail

Collisions::Collisions(
        DDomSpTor3DV2D dom_sp_tor3D_v2D,
        device_t<DViewMu> coeff_intdmu,
        device_t<DViewVpar> coeff_intdvpar,
        device_t<DViewTor1> nustar0_r)
    : m_operator_handle {}
    , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
    , m_hat_As {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_As"), ddc::select<Species>(dom_sp_tor3D_v2D).size()}   // m_hat_As{"m_hat_As",ddc::select<Species>(dom_sp_tor3D_v2D)}
    , m_hat_Zs {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_Zs"), ddc::select<Species>(dom_sp_tor3D_v2D).size()}
    , m_rg {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_rg"), ddc::select<GridTor1>(dom_sp_tor3D_v2D).size()}
    , m_q_r {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_q_r"), ddc::select<GridTor1>(dom_sp_tor3D_v2D).size()}
    , m_mask_buffer_r {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_buffer_r"), ddc::select<GridTor1>(dom_sp_tor3D_v2D).size()}
    , m_mask_LIM {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_LIM"), ddc::select<GridTor1>(dom_sp_tor3D_v2D).size(), ddc::select<GridTor2>(dom_sp_tor3D_v2D).size()}
    , m_B_norm {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_B_norm"), ddc::select<GridTor1>(dom_sp_tor3D_v2D).size(), ddc::select<GridTor2>(dom_sp_tor3D_v2D).size()}
    , m_Bstar_s {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_Bstar_s"),
              ddc::select<GridVpar>(dom_sp_tor3D_v2D).size(),
              ddc::select<GridTor1>(dom_sp_tor3D_v2D).size(),
              ddc::select<GridTor2>(dom_sp_tor3D_v2D).size(),
              ddc::select<Species>(dom_sp_tor3D_v2D).size()}
    , m_mug {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mug"),
              ddc::select<GridMu>(dom_sp_tor3D_v2D).size()}
    , m_vparg {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_vparg"),
              ddc::select<GridVpar>(dom_sp_tor3D_v2D).size()}
{
    detail::DoCombMatComputation(m_comb_mat);

    Kokkos::deep_copy(m_hat_As, 1.0); //ddc::fill()
    Kokkos::deep_copy(m_hat_Zs, 1.0);
    Kokkos::deep_copy(m_rg, 1.0);
    Kokkos::deep_copy(m_q_r, 1.0);
    Kokkos::deep_copy(m_mask_buffer_r, 0.0); // Masked if >= 0.99
    Kokkos::deep_copy(m_mask_LIM, 0.0); // Masked if >= 0.99
    Kokkos::deep_copy(m_B_norm, 1.0);
    Kokkos::deep_copy(m_Bstar_s, 1.0);

    detail::DumpCoordinates(m_mug, ddc::select<GridMu>(dom_sp_tor3D_v2D));
    detail::DumpCoordinates(m_vparg, ddc::select<GridVpar>(dom_sp_tor3D_v2D));

    // NOTE: We need to fence because DumpCoordinates is asynchronous on GPUs.
    Kokkos::fence();

    m_operator_handle = detail::DoOperatorInitialization(
            ddc::select<GridMu>(dom_sp_tor3D_v2D).size(),
            ddc::select<GridVpar>(dom_sp_tor3D_v2D).size(),
            ddc::select<GridTor1>(dom_sp_tor3D_v2D).size(),
            ddc::select<GridTor2>(dom_sp_tor3D_v2D).size(),
            ddc::select<GridTor3>(dom_sp_tor3D_v2D).size(),
            ddc::select<Species>(dom_sp_tor3D_v2D).size(),
            /* the_local_domain_r_offset */ 0 + ddc::select<GridTor1>(dom_sp_tor3D_v2D).size() - 1,
            m_mug.data(),
            m_vparg.data(),
            coeff_intdmu.data_handle(),
            coeff_intdvpar.data_handle(),
            nustar0_r.data_handle(),
            m_comb_mat.data(),
            m_hat_As.data(), // m_hat_As.data_handle()   --> pointer
            m_hat_Zs.data(),
            m_rg.data(),
            m_q_r.data(),
            m_mask_buffer_r.data(),
            m_mask_LIM.data(),
            m_B_norm.data(),
            m_Bstar_s.data());
}

Collisions::~Collisions()
{
    detail::DoOperatorDeinitialization(static_cast<::koliop_Operator>(m_operator_handle));
}

void Collisions::operator()(device_t<DSpanSpTor3DV2D> all_f_distribution, double deltat_coll) const
{
    if (::koliop_Collision(
                static_cast<::koliop_Operator>(m_operator_handle),
                deltat_coll,
                all_f_distribution.data_handle())
        != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }

    // NOTE: While gslx is not fully stream compatible, just fence at the end of
    // the computation.
    if (::koliop_Fence(static_cast<::koliop_Operator>(m_operator_handle))
        != KOLIOP_STATUS_SUCCESS) {
        GSLX_ASSERT(false);
    }
}
