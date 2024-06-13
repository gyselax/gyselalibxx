#pragma once

#include <ddc/ddc.hpp>

#include <assert.hpp>
#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "koliop_interface.hpp"

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
template <class FDistribDomain, class GridVpar, class GridMu>
class CollisionSpVparMu /* : public IRightHandSide */
{
private:
    using Species = IDimSp;
    using fdistrib_domain_tags = ddc::to_type_seq_t<FDistribDomain>;
    // Ensure expected types appear in distribution domain
    static_assert(
            ddc::in_tags_v<Species, fdistrib_domain_tags>,
            "Species is missing from distribution function domain");
    static_assert(
            ddc::in_tags_v<GridVpar, fdistrib_domain_tags>,
            "Vpar is missing from distribution function domain");
    static_assert(
            ddc::in_tags_v<GridMu, fdistrib_domain_tags>,
            "Mu is missing from distribution function domain");
    // Ensure expected types appear in distribution domain at expected position
    static_assert(
            ddc::type_seq_rank_v<Species, fdistrib_domain_tags> == 0,
            "Species should appear first in the distribution function domain");
    static_assert(
            ddc::type_seq_rank_v<GridMu, fdistrib_domain_tags> == (FDistribDomain::rank() - 1),
            "Mu should appear last in the distribution function domain");
    static_assert(
            ddc::type_seq_rank_v<GridVpar, fdistrib_domain_tags> == (FDistribDomain::rank() - 2),
            "Vpar should appear second to last in the distribution function domain");

public:
    /// Type alias for the domain of the magnetic moment.
    using DDomMu = ddc::DiscreteDomain<GridMu>;
    /// Type alias for the domain of the velocity parallel to the magnetic field.
    using DDomVpar = ddc::DiscreteDomain<GridVpar>;
    /// Type alias for a constant reference to a Chunk on GPU defined on a grid of magnetic moments.
    using DViewMu = ddc::ChunkSpan<
            double const,
            DDomMu,
            std::experimental::layout_right,
            Kokkos::DefaultExecutionSpace::
                    memory_space>; // Equivalent to device_t<ddc::ChunkSpan<double const, DDomMu>>
    /// Type alias for a constant reference to a Chunk on GPU defined on a grid of parallel velocities.
    using DViewVpar = ddc::ChunkSpan<
            double const,
            DDomVpar,
            std::experimental::layout_right,
            Kokkos::DefaultExecutionSpace::
                    memory_space>; // Equivalent to device_t<ddc::ChunkSpan<double const, DDomMu>>

    /// Type alias for a reference to the distribution function stored on GPU.
    using FDistribSpan = device_t<ddc::ChunkSpan<double, FDistribDomain>>;

private:
    // Shortcut to get the type at a specified index in the Chunk domain
    template <std::size_t I, class ChunkType>
    using domain_element_t
            = ddc::type_seq_element_t<I, ddc::to_type_seq_t<typename ChunkType::mdomain_type>>;

    // Get the size of a chunk. This function returns 1 if a double is passed.
    template <class ChunkType>
    static std::size_t get_chunk_size(ChunkType chunk_1d)
    {
        static_assert(std::is_same_v<ChunkType, double> || ddc::is_borrowed_chunk_v<ChunkType>);
        if constexpr (ddc::is_chunk_v<ChunkType>) {
            return chunk_1d.domain().size();
        } else {
            return 1;
        }
    }

    // Get the size of the radial domain from nustar0_r
    template <class NuStar0Type>
    static std::size_t get_r_domain_size(NuStar0Type nustar0_r)
    {
        static_assert(
                std::is_same_v<
                        NuStar0Type,
                        double> || (ddc::is_borrowed_chunk_v<NuStar0Type> && NuStar0Type::rank() == 1));
        return get_chunk_size(nustar0_r);
    }

    // Get the size of the poloidal domain from B_norm by excluding the size of nustar0_r
    template <class NuStar0Type, class BNormType>
    static std::size_t get_theta_domain_size(NuStar0Type nustar0_r, BNormType B_norm)
    {
        std::size_t const r_size = get_chunk_size(nustar0_r);
        std::size_t const r_theta_size = get_chunk_size(B_norm);
        return r_theta_size / r_size;
    }

public:
    /**
     * @brief Create instance of CollisionSpVparMu class
     *
     * @param[in] fdistrib_domain
     *      domain (species, 3D space, 2D velocity space) on which the collision operator acts
     *      ATTENTION it must contain the all domain in species and 2D velocity.
     * @param[in] coeff_intdmu
     *      quadrature coefficient in mu direction
     * @param[in] coeff_intdvpar
     *      quadrature coefficient in vpar direction
     * @param[in] nustar0_r
     *      radial profile of collisionality
     * @param[in] B_norm
     *      magnetic field norm in (r,theta)
     */
    template <class NuStar0Type, class BNormType>
    CollisionSpVparMu(
            FDistribDomain fdistrib_domain,
            DViewMu coeff_intdmu,
            DViewVpar coeff_intdvpar,
            NuStar0Type nustar0_r,
            BNormType B_norm)
    : m_operator_handle {}
    , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
    , m_hat_As {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_As"), ddc::select<Species>(fdistrib_domain).size()}   // m_hat_As{"m_hat_As",ddc::select<Species>(fdistrib_domain)}
    , m_hat_Zs {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_Zs"), ddc::select<Species>(fdistrib_domain).size()}
    , m_rg {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_rg"), get_r_domain_size(nustar0_r)}
    , m_q_r {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_q_r"), get_r_domain_size(nustar0_r)}
    , m_mask_buffer_r {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_buffer_r"), get_r_domain_size(nustar0_r)}
    , m_mask_LIM {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_LIM"), get_r_domain_size(nustar0_r), get_theta_domain_size(nustar0_r, B_norm)}
    , m_B_norm {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_B_norm"), get_r_domain_size(nustar0_r), get_theta_domain_size(nustar0_r, B_norm)}
    , m_Bstar_s {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_Bstar_s"),
              ddc::select<GridVpar>(fdistrib_domain).size(),
              get_r_domain_size(nustar0_r),
              get_theta_domain_size(nustar0_r, B_norm),
              ddc::select<Species>(fdistrib_domain).size()}
    , m_mug {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mug"),
              ddc::select<GridMu>(fdistrib_domain).size()}
    , m_vparg {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_vparg"),
              ddc::select<GridVpar>(fdistrib_domain).size()}
    {
        // Check that the distribution function is correctly ordered
        if constexpr (ddc::is_chunk_v<NuStar0Type>) {
            using GridR = domain_element_t<0, NuStar0Type>;
            if constexpr (ddc::is_chunk_v<BNormType>) {
                static_assert(std::is_same_v<GridR, domain_element_t<1, BNormType>>);
                using GridTheta = domain_element_t<0, BNormType>;
                static_assert(
                        ddc::type_seq_rank_v<
                                GridR,
                                fdistrib_domain_tags> == (FDistribDomain::rank() - 3),
                        "R should appear third to last in the distribution function domain");
                static_assert(
                        ddc::type_seq_rank_v<
                                GridTheta,
                                fdistrib_domain_tags> == (FDistribDomain::rank() - 4),
                        "Theta should appear fourth to last in the distribution function domain");
            } else {
                static_assert(
                        ddc::type_seq_rank_v<
                                GridR,
                                fdistrib_domain_tags> == (FDistribDomain::rank() - 3),
                        "R should appear third to last in the distribution function domain");
            }
        } else if constexpr (ddc::is_chunk_v<BNormType>) {
            using GridTheta = domain_element_t<0, BNormType>;
            static_assert(
                    ddc::type_seq_rank_v<
                            GridTheta,
                            fdistrib_domain_tags> == (FDistribDomain::rank() - 3),
                    "Theta should appear third to last in the distribution function domain");
        }
        koliop_interface::DoCombMatComputation(m_comb_mat);

        Kokkos::deep_copy(m_hat_As, 1.0); //ddc::fill()
        Kokkos::deep_copy(m_hat_Zs, 1.0);
        Kokkos::deep_copy(m_rg, 1.0);
        Kokkos::deep_copy(m_q_r, 1.0);
        Kokkos::deep_copy(m_mask_buffer_r, 0.0); // Masked if >= 0.99
        Kokkos::deep_copy(m_mask_LIM, 0.0); // Masked if >= 0.99
        Kokkos::deep_copy(m_B_norm, 1.0); // To replace and use argument instead
        Kokkos::deep_copy(m_Bstar_s, 1.0);

        koliop_interface::DumpCoordinates(m_mug, ddc::select<GridMu>(fdistrib_domain));
        koliop_interface::DumpCoordinates(m_vparg, ddc::select<GridVpar>(fdistrib_domain));

        // NOTE: We need to fence because DumpCoordinates is asynchronous on GPUs.
        Kokkos::fence();

        std::size_t const n_mu = ddc::select<GridMu>(fdistrib_domain).size();
        std::size_t const n_vpar = ddc::select<GridVpar>(fdistrib_domain).size();
        std::size_t const n_r = get_r_domain_size(nustar0_r);
        std::size_t const n_theta = get_theta_domain_size(nustar0_r, B_norm);
        std::size_t const n_sp = ddc::select<Species>(fdistrib_domain).size();
        std::size_t const n_batch = fdistrib_domain.size() / (n_mu * n_vpar * n_r * n_theta * n_sp);

        double* nustar0_r_ptr;
        if constexpr (std::is_same_v<NuStar0Type, double>) {
            nustar0_r_ptr = &nustar0_r;
        } else {
            nustar0_r_ptr = nustar0_r.data_handle();
        }

        m_operator_handle = koliop_interface::DoOperatorInitialization(
                n_mu,
                n_vpar,
                n_r,
                n_theta,
                n_batch,
                n_sp,
                /* the_local_domain_r_offset */ 0 + n_r - 1,
                m_mug.data(),
                m_vparg.data(),
                coeff_intdmu.data_handle(),
                coeff_intdvpar.data_handle(),
                nustar0_r_ptr,
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

    ~CollisionSpVparMu()
    {
        koliop_interface::DoOperatorDeinitialization(
                static_cast<::koliop_Operator>(m_operator_handle));
    }

    /**
     * @brief Apply the collision operator to the distribution functions of all species on all species
     *
     * @param[inout] all_f_distribution
     *      All distribution functions
     * @param[in] deltat_coll
     *      Collision time step
     */
    void operator()(FDistribSpan all_f_distribution, double deltat_coll) const
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

protected:
    /**
     * Opaque type representing the operator (due to the C interface)
    */
    ::koliop_Operator m_operator_handle;

    // NOTE: Some of these arrays should come from a parent class or manager.
    // They resides in this class while we wait for their implementation.
    /// Combinatory (6x6) matrix computed only one times at initialisation. Rk: 6 = 2*(Npolmax-1) + 1 + 1
    koliop_interface::MDL<double[6][6]> m_comb_mat;
    /// Normalized masses for all species
    koliop_interface::MDL<double*> m_hat_As; // DFieldSp
    /// Normalized charges for all species
    koliop_interface::MDL<double*> m_hat_Zs;
    /// Mesh points in the radial direction
    // TODO: See if we need m_rg ?
    koliop_interface::MDL<double*> m_rg;
    /// Radial safety factor profile
    koliop_interface::MDL<double*> m_q_r;
    /// Mask used to avoid to apply collision in certain region
    // TODO: Understand this mask
    koliop_interface::MDL<double*> m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    koliop_interface::MDL<double**> m_mask_LIM;
    /// B norm in (r,theta)
    // TODO: Attention this must be 3D --> transfer it in a 1D array ?
    koliop_interface::MDL<double**> m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // TODO : Must be 5D for full 3D geometry
    koliop_interface::MDL<double****> m_Bstar_s;
    /// grid in mu direction
    koliop_interface::MDL<double*> m_mug;
    /// grid in vpar direction
    koliop_interface::MDL<double*> m_vparg;
};
