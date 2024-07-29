#pragma once

#include <ddc/ddc.hpp>

#include <assert.hpp>
#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "collisions_dimensions.hpp"
#include "koliop_interface.hpp"

/**
 * @brief A class which computes the collision operator in (vpar,mu)
 */
template <
        class CollInfo,
        class FDistribDomain,
        class GridVpar,
        class GridMu,
        class InputDFieldThetaR>
class CollisionSpVparMu /* : public IRightHandSide */
{
private:
    using InputDFieldR = typename CollInfo::radial_chunk_type;
    using fdistrib_domain_tags = ddc::to_type_seq_t<FDistribDomain>;
    using GridR = typename collisions_dimensions::ExtractRDim<InputDFieldR>::type;
    using GridTheta =
            typename collisions_dimensions::ExtractThetaDim<InputDFieldThetaR, GridR>::type;

    // Validate template types
    static_assert(ddc::is_discrete_domain_v<FDistribDomain>);
    static_assert(FDistribDomain::rank() >= 3 && FDistribDomain::rank() <= 6);
    static_assert((std::is_same_v<InputDFieldR, double>) || ddc::is_borrowed_chunk_v<InputDFieldR>);
    static_assert(
            (std::is_same_v<InputDFieldThetaR, double>)
            || (ddc::is_borrowed_chunk_v<InputDFieldThetaR>));
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
    // Ensure that uniform
    // [TODO] Restore it as soon as the geometry5D is deleted
    //static_assert(ddc::is_uniform_point_sampling_v<GridVpar>);

    static_assert(
            ddc::type_seq_rank_v<Species, fdistrib_domain_tags> == 0,
            "Species should appear first in the distribution function domain");
    static_assert(
            collisions_dimensions::
                    order_of_last_grids<fdistrib_domain_tags, GridTheta, GridR, GridVpar, GridMu>(),
            "Misordered distribution function domain detected. Koliop expects (Sp, Phi, Theta, R, "
            "Vpar, Mu)");

public:
    /// Type alias for the domain of the radial points
    using DDomR = ddc::DiscreteDomain<GridR>;
    /// Type alias for the domain of the poloidal plane
    using DDomThetaR = ddc::DiscreteDomain<GridTheta, GridR>;
    /// Type alias for the domain containing (species, theta, r, vpar)
    using DDomSpThetaRVpar = ddc::DiscreteDomain<Species, GridTheta, GridR, GridVpar>;

    /// Type alias for the domain of the magnetic moment.
    using DDomMu = ddc::DiscreteDomain<GridMu>;
    /// Type alias for the domain of the velocity parallel to the magnetic field.
    using DDomVpar = ddc::DiscreteDomain<GridVpar>;
    /// Type alias for a field on a grid of species
    using DFieldMemSp = device_t<ddc::Chunk<double, IdxRangeSp>>;
    /// Type alias for a field on a grid of radial values
    using DFieldR = device_t<ddc::Chunk<double, DDomR>>;
    /// Type alias for a field on a grid of magnetic moments
    using DFieldMu = device_t<ddc::Chunk<double, DDomMu>>;
    /// Type alias for a field on a grid of parallel velocities
    using DFieldVpar = device_t<ddc::Chunk<double, DDomVpar>>;
    /// Type alias for a field on a grid on a poloidal plane
    using DFieldThetaR = device_t<ddc::Chunk<double, DDomThetaR>>;
    /// Type alias for a field on a grid of species, poloidal plane and parallel velocities
    using DFieldSpThetaRVpar = device_t<ddc::Chunk<double, DDomSpThetaRVpar>>;
    /// Type alias for a span of a field defined on a grid of radial values
    using DSpanR = device_t<ddc::ChunkSpan<double, DDomR>>;
    /// Type alias for a span of a field defined on a grid on a poloidal plane
    using DSpanThetaR = device_t<ddc::ChunkSpan<double, DDomThetaR>>;
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
    /**
     * @brief Copy the information in src (received as an input to the class) into a 1D radial profile.
     *
     * @param[in] src The source from which the data is copied.
     * @param[out] dst The destination into which the data is copied.
     */
    void deepcopy_radial_profile(DSpanR dst, InputDFieldR src)
    {
        if constexpr (collisions_dimensions::is_spoofed_dim_v<GridR>) {
            ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), dst, src);
        } else {
            ddc::parallel_deepcopy(dst, src);
        }
    }

public:
    /**
     * @brief Copy the information in src (received as an input to the class) into a 2D field defined
     * over (r, theta).
     * This function should be private but is public due to Kokkos restrictions.
     *
     * @param[in] src The source from which the data is copied.
     * @param[out] dst The destination into which the data is copied.
     */
    void deepcopy_poloidal_plane(DSpanThetaR dst, InputDFieldThetaR src)
    {
        if constexpr ((collisions_dimensions::is_spoofed_dim_v<GridR>)&&(
                              collisions_dimensions::is_spoofed_dim_v<GridTheta>)) {
            ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), dst, src);
        } else if constexpr ((!collisions_dimensions::is_spoofed_dim_v<GridR>)&&(
                                     !collisions_dimensions::is_spoofed_dim_v<GridTheta>)) {
            ddc::parallel_deepcopy(dst, src);
        } else {
            using NonSpoofDim = std::
                    conditional_t<collisions_dimensions::is_spoofed_dim_v<GridR>, GridTheta, GridR>;
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    dst.domain(),
                    KOKKOS_LAMBDA(ddc::DiscreteElement<GridTheta, GridR> idx) {
                        dst(idx) = src(ddc::select<NonSpoofDim>(idx));
                    });
        }
    }

public:
    /**
     * @brief Create instance of CollisionSpVparMu class
     *
     * @param[in] collision_info
     *      class containing rg, safety_factor and nustar0
     * @param[in] fdistrib_domain
     *      domain (species, 3D space, 2D velocity space) on which the collision operator acts
     *      ATTENTION it must contain the all domain in species and 2D velocity.
     * @param[in] coeff_intdmu
     *      quadrature coefficient in mu direction
     * @param[in] coeff_intdvpar
     *      quadrature coefficient in vpar direction
     * @param[in] B_norm
     *      magnetic field norm in (r,theta)
     */
    CollisionSpVparMu(
            CollInfo const& collision_info,
            FDistribDomain fdistrib_domain,
            DViewMu coeff_intdmu,
            DViewVpar coeff_intdvpar,
            InputDFieldThetaR B_norm)
        : m_operator_handle {}
        , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
        , m_hat_As {"m_hat_As", collisions_dimensions::get_idx_range<Species>(fdistrib_domain)}
        , m_hat_Zs {"m_hat_Zs", collisions_dimensions::get_idx_range<Species>(fdistrib_domain)}
        , m_nustar0_r {"m_nustar0_r", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_rg {"m_rg", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_safety_factor {"m_safety_factor", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_mask_buffer_r {"m_mask_buffer_r", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_mask_LIM {"m_mask_LIM", DDomThetaR {collisions_dimensions::get_idx_range<GridTheta, GridR>(fdistrib_domain)}}
        , m_B_norm {"m_B_norm", DDomThetaR {collisions_dimensions::get_idx_range<GridTheta, GridR>(fdistrib_domain)}}
        , m_Bstar_s {"m_Bstar_s", DDomSpThetaRVpar {collisions_dimensions::get_idx_range<Species, GridTheta, GridR, GridVpar>(fdistrib_domain)}}
        , m_mug {"m_mug", ddc::select<GridMu>(fdistrib_domain)}
        , m_vparg {"m_vparg", ddc::select<GridVpar>(fdistrib_domain)}
    {
        // Check that the distribution function is correctly ordered
        koliop_interface::DoCombMatComputation(m_comb_mat);

        IdxRangeSp idxrange_sp = ddc::select<Species>(fdistrib_domain);
        // --> Initialize the mass species
        ddc::ChunkSpan hat_As_host = ddc::discrete_space<Species>().masses()[idxrange_sp];
        ddc::parallel_deepcopy(m_hat_As.span_view(), hat_As_host);
        // --> Initialize the charge species
        ddc::ChunkSpan hat_Zs_host = ddc::discrete_space<Species>().charges()[idxrange_sp];
        ddc::parallel_deepcopy(m_hat_Zs.span_view(), hat_Zs_host);

        // --> Initialize the other quantities needed in koliop
        // TODO: Put Bstar_s as an input variable of the constructor (something more specific than what is done for B_norm must be done)
        ddc::parallel_fill(m_Bstar_s.span_view(), 1.0);

        // --> Initialization of the masks that have no sense here to 0.
        ddc::parallel_fill(m_mask_buffer_r.span_view(), 0.0); // Masked if >= 0.99
        ddc::parallel_fill(m_mask_LIM.span_view(), 0.0); // Masked if >= 0.99

        // --> Initialization of vpar and mu grids
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), m_mug.span_view());
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), m_vparg.span_view());

        // NOTE: We need to fence because DumpCoordinates is asynchronous on GPUs.
        Kokkos::fence();

        std::size_t const n_mu = ddc::select<GridMu>(fdistrib_domain).size();
        std::size_t const n_vpar = ddc::select<GridVpar>(fdistrib_domain).size();
        std::size_t const n_r = collisions_dimensions::get_idx_range<GridR>(fdistrib_domain).size();
        std::size_t const n_theta
                = collisions_dimensions::get_idx_range<GridTheta>(fdistrib_domain).size();
        std::size_t const n_sp = ddc::select<Species>(fdistrib_domain).size();
        std::size_t const n_batch = fdistrib_domain.size() / (n_mu * n_vpar * n_r * n_theta * n_sp);

        deepcopy_radial_profile(m_rg.span_view(), collision_info.rg());
        deepcopy_radial_profile(m_safety_factor.span_view(), collision_info.safety_factor());
        deepcopy_radial_profile(m_nustar0_r.span_view(), collision_info.nustar0());
        deepcopy_poloidal_plane(m_B_norm.span_view(), B_norm);

        m_operator_handle = koliop_interface::DoOperatorInitialization(
                n_mu,
                n_vpar,
                n_r,
                n_theta,
                n_batch,
                n_sp,
                collision_info.collisions_interspecies(),
                /* the_local_domain_r_offset */ 0 + n_r - 1,
                m_mug.data_handle(),
                m_vparg.data_handle(),
                coeff_intdmu.data_handle(),
                coeff_intdvpar.data_handle(),
                m_nustar0_r.data_handle(),
                m_comb_mat.data(),
                m_hat_As.data_handle(),
                m_hat_Zs.data_handle(),
                m_rg.data_handle(),
                m_safety_factor.data_handle(),
                m_mask_buffer_r.data_handle(),
                m_mask_LIM.data_handle(),
                m_B_norm.data_handle(),
                m_Bstar_s.data_handle());
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
    DFieldMemSp m_hat_As;
    /// Normalized charges for all species
    DFieldMemSp m_hat_Zs;
    /// Radial profile of nustar0_r
    DFieldR m_nustar0_r;
    /// Mesh points in the radial direction
    // [TODO]: See if we need m_rg ?
    DFieldR m_rg;
    /// Radial safety factor profile
    DFieldR m_safety_factor;
    /// Mask used to avoid to apply collision in certain region
    // [TODO]: This mask should maybe be deleted in C++ version
    DFieldR m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    DFieldThetaR m_mask_LIM;
    /// B norm in (r,theta)
    // [TODO] Attention this must be 3D for generalization to 3D geometry--> transfer it in a 1D array ?
    DFieldThetaR m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // [TODO] Must be 5D for full 3D geometry
    DFieldSpThetaRVpar m_Bstar_s;
    /// grid in mu direction
    DFieldMu m_mug;
    /// grid in vpar direction
    DFieldVpar m_vparg;
};
