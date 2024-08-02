#pragma once
#include <ddc/ddc.hpp>

#include <assert.hpp>
#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "collisions_dimensions.hpp"
#include "ddc_aliases.hpp"
#include "koliop_interface.hpp"

/**
 * @brief A class which computes the collision operator in (vpar,mu)
 */
template <
        class CollInfo,
        class FDistribIdxRange,
        class GridVpar,
        class GridMu,
        class InputDFieldThetaR>
class CollisionSpVparMu /* : public IRightHandSide */
{
private:
    using InputDFieldR = typename CollInfo::radial_chunk_type;
    using fdistrib_grids = ddc::to_type_seq_t<FDistribIdxRange>;
    using GridR = typename collisions_dimensions::ExtractRDim<InputDFieldR>::type;
    using GridTheta =
            typename collisions_dimensions::ExtractThetaDim<InputDFieldThetaR, GridR>::type;

    // Validate template types
    static_assert(ddc::is_discrete_domain_v<FDistribIdxRange>);
    static_assert(FDistribIdxRange::rank() >= 3 && FDistribIdxRange::rank() <= 6);
    static_assert((std::is_same_v<InputDFieldR, double>) || ddc::is_borrowed_chunk_v<InputDFieldR>);
    static_assert(
            (std::is_same_v<InputDFieldThetaR, double>)
            || (ddc::is_borrowed_chunk_v<InputDFieldThetaR>));
    // Ensure expected types appear in distribution index range
    static_assert(
            ddc::in_tags_v<Species, fdistrib_grids>,
            "The distribution function must be defined over the species");
    static_assert(
            ddc::in_tags_v<GridVpar, fdistrib_grids>,
            "The distribution function must be defined over the Vpar direction");
    static_assert(
            ddc::in_tags_v<GridMu, fdistrib_grids>,
            "The distribution function must be defined over the Mu direction");
    // Ensure that uniform
    // [TODO] Restore it as soon as the geometry5D is deleted
    //static_assert(ddc::is_uniform_point_sampling_v<GridVpar>);

    static_assert(
            ddc::type_seq_rank_v<Species, fdistrib_grids> == 0,
            "Species should be the first index in the distribution function index set");
    static_assert(
            collisions_dimensions::
                    order_of_last_grids<fdistrib_grids, GridTheta, GridR, GridVpar, GridMu>(),
            "Misordered inex set for the distribution function. Koliop expects (Sp, Phi, Theta, R, "
            "Vpar, Mu)");

public:
    /// Type alias for the index range of the radial points
    using IdxRangeR = IdxRange<GridR>;
    /// Type alias for the index range of the poloidal plane
    using IdxRangeThetaR = IdxRange<GridTheta, GridR>;
    /// Type alias for the index range containing (species, theta, r, vpar)
    using IdxRangeSpThetaRVpar = IdxRange<Species, GridTheta, GridR, GridVpar>;

    /// Type alias for the index range of the magnetic moment.
    using IdxRangeMu = IdxRange<GridMu>;
    /// Type alias for the index range of the velocity parallel to the magnetic field.
    using IdxRangeVpar = IdxRange<GridVpar>;
    /// Type alias for a field on a grid of species
    using DFieldMemSp = FieldMem<double, IdxRangeSp>;
    /// Type alias for a field on a grid of radial values
    using DFieldMemR = FieldMem<double, IdxRangeR>;
    /// Type alias for a field on a grid of magnetic moments
    using DFieldMemMu = FieldMem<double, IdxRangeMu>;
    /// Type alias for a field on a grid of parallel velocities
    using DFieldMemVpar = FieldMem<double, IdxRangeVpar>;
    /// Type alias for a field on a grid on a poloidal plane
    using DFieldMemThetaR = FieldMem<double, IdxRangeThetaR>;
    /// Type alias for a field on a grid of species, poloidal plane and parallel velocities
    using DFieldMemSpThetaRVpar = FieldMem<double, IdxRangeSpThetaRVpar>;
    /// Type alias for a span of a field defined on a grid of radial values
    using DFieldR = Field<double, IdxRangeR>;
    /// Type alias for a span of a field defined on a grid on a poloidal plane
    using DFieldThetaR = Field<double, IdxRangeThetaR>;
    /// Type alias for a constant reference to a Chunk on GPU defined on a grid of magnetic moments.
    using DConstFieldMu
            = Field<double const,
                    IdxRangeMu,
                    std::experimental::layout_right,
                    Kokkos::DefaultExecutionSpace::
                            memory_space>; // Equivalent to Field<double const, IdxRangeMu>
    /// Type alias for a constant reference to a Chunk on GPU defined on a grid of parallel velocities.
    using DConstFieldVpar
            = Field<double const,
                    IdxRangeVpar,
                    std::experimental::layout_right,
                    Kokkos::DefaultExecutionSpace::
                            memory_space>; // Equivalent to Field<double const, IdxRangeMu>

    /// Type alias for a reference to the distribution function stored on GPU.
    using FDistribField = Field<double, FDistribIdxRange>;

private:
    /**
     * @brief Copy the information in src (received as an input to the class) into a 1D radial profile.
     *
     * @param[in] src The source from which the data is copied.
     * @param[out] dst The destination into which the data is copied.
     */
    void deepcopy_radial_profile(DFieldR dst, InputDFieldR src)
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
    void deepcopy_poloidal_plane(DFieldThetaR dst, InputDFieldThetaR src)
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
                    get_idx_range(dst),
                    KOKKOS_LAMBDA(Idx<GridTheta, GridR> idx) {
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
     * @param[in] fdistrib_idx_range
     *      index range (species, 3D space, 2D velocity space) on which the collision operator acts
     *      ATTENTION it must contain the all index range in species and 2D velocity.
     * @param[in] coeff_intdmu
     *      quadrature coefficient in mu direction
     * @param[in] coeff_intdvpar
     *      quadrature coefficient in vpar direction
     * @param[in] B_norm
     *      magnetic field norm in (r,theta)
     */
    CollisionSpVparMu(
            CollInfo const& collision_info,
            FDistribIdxRange fdistrib_idx_range,
            DConstFieldMu coeff_intdmu,
            DConstFieldVpar coeff_intdvpar,
            InputDFieldThetaR B_norm)
        : m_operator_handle {}
        , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
        , m_hat_As {"m_hat_As", collisions_dimensions::get_idx_range<Species>(fdistrib_idx_range)}
        , m_hat_Zs {"m_hat_Zs", collisions_dimensions::get_idx_range<Species>(fdistrib_idx_range)}
        , m_nustar0_r {"m_nustar0_r", collisions_dimensions::get_idx_range<GridR>(fdistrib_idx_range)}
        , m_rg {"m_rg", collisions_dimensions::get_idx_range<GridR>(fdistrib_idx_range)}
        , m_safety_factor {"m_safety_factor", collisions_dimensions::get_idx_range<GridR>(fdistrib_idx_range)}
        , m_mask_buffer_r {"m_mask_buffer_r", collisions_dimensions::get_idx_range<GridR>(fdistrib_idx_range)}
        , m_mask_LIM {"m_mask_LIM", IdxRangeThetaR {collisions_dimensions::get_idx_range<GridTheta, GridR>(fdistrib_idx_range)}}
        , m_B_norm {"m_B_norm", IdxRangeThetaR {collisions_dimensions::get_idx_range<GridTheta, GridR>(fdistrib_idx_range)}}
        , m_Bstar_s {"m_Bstar_s", IdxRangeSpThetaRVpar {collisions_dimensions::get_idx_range<Species, GridTheta, GridR, GridVpar>(fdistrib_idx_range)}}
        , m_mug {"m_mug", ddc::select<GridMu>(fdistrib_idx_range)}
        , m_vparg {"m_vparg", ddc::select<GridVpar>(fdistrib_idx_range)}
    {
        // Check that the distribution function is correctly ordered
        koliop_interface::DoCombMatComputation(m_comb_mat);

        IdxRangeSp idxrange_sp = ddc::select<Species>(fdistrib_idx_range);
        // --> Initialize the mass species
        ddc::ChunkSpan hat_As_host = ddc::discrete_space<Species>().masses()[idxrange_sp];
        ddc::parallel_deepcopy(get_field(m_hat_As), hat_As_host);
        // --> Initialize the charge species
        ddc::ChunkSpan hat_Zs_host = ddc::discrete_space<Species>().charges()[idxrange_sp];
        ddc::parallel_deepcopy(get_field(m_hat_Zs), hat_Zs_host);

        // --> Initialize the other quantities needed in koliop
        // TODO: Put Bstar_s as an input variable of the constructor (something more specific than what is done for B_norm must be done)
        ddc::parallel_fill(get_field(m_Bstar_s), 1.0);

        // --> Initialization of the masks that have no sense here to 0.
        ddc::parallel_fill(get_field(m_mask_buffer_r), 0.0); // Masked if >= 0.99
        ddc::parallel_fill(get_field(m_mask_LIM), 0.0); // Masked if >= 0.99

        // --> Initialization of vpar and mu grids
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_mug));
        ddcHelper::dump_coordinates(Kokkos::DefaultExecutionSpace(), get_field(m_vparg));

        // NOTE: We need to fence because DumpCoordinates is asynchronous on GPUs.
        Kokkos::fence();

        std::size_t const n_mu = ddc::select<GridMu>(fdistrib_idx_range).size();
        std::size_t const n_vpar = ddc::select<GridVpar>(fdistrib_idx_range).size();
        std::size_t const n_r
                = collisions_dimensions::get_idx_range<GridR>(fdistrib_idx_range).size();
        std::size_t const n_theta
                = collisions_dimensions::get_idx_range<GridTheta>(fdistrib_idx_range).size();
        std::size_t const n_sp = ddc::select<Species>(fdistrib_idx_range).size();
        std::size_t const n_batch
                = fdistrib_idx_range.size() / (n_mu * n_vpar * n_r * n_theta * n_sp);

        deepcopy_radial_profile(get_field(m_rg), collision_info.rg());
        deepcopy_radial_profile(get_field(m_safety_factor), collision_info.safety_factor());
        deepcopy_radial_profile(get_field(m_nustar0_r), collision_info.nustar0());
        deepcopy_poloidal_plane(get_field(m_B_norm), B_norm);

        m_operator_handle = koliop_interface::DoOperatorInitialization(
                n_mu,
                n_vpar,
                n_r,
                n_theta,
                n_batch,
                n_sp,
                collision_info.collisions_interspecies(),
                /* the_local_index range_r_offset */ 0 + n_r - 1,
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
    void operator()(FDistribField all_f_distribution, double deltat_coll) const
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
    DFieldMemR m_nustar0_r;
    /// Mesh points in the radial direction
    // [TODO]: See if we need m_rg ?
    DFieldMemR m_rg;
    /// Radial safety factor profile
    DFieldMemR m_safety_factor;
    /// Mask used to avoid to apply collision in certain region
    // [TODO]: This mask should maybe be deleted in C++ version
    DFieldMemR m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    DFieldMemThetaR m_mask_LIM;
    /// B norm in (r,theta)
    // [TODO] Attention this must be 3D for generalization to 3D geometry--> transfer it in a 1D array ?
    DFieldMemThetaR m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // [TODO] Must be 5D for full 3D geometry
    DFieldMemSpThetaRVpar m_Bstar_s;
    /// grid in mu direction
    DFieldMemMu m_mug;
    /// grid in vpar direction
    DFieldMemVpar m_vparg;
};
