#pragma once

#include <ddc/ddc.hpp>

#include <assert.hpp>
#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "koliop_interface.hpp"

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
    // Ensure that uniform
    // [TODO] Restore it as soon as the geometry5D is deleted
    //static_assert(ddc::is_uniform_point_sampling_v<GridVpar>);

public:
    /// Type alias for a field on species domain
    using DFieldSp = device_t<ddc::Chunk<double, IDomainSp>>;
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

    // Get the size of the radial domain from nustar
    template <class DFieldR>
    static std::size_t get_r_domain_size(DFieldR nustar)
    {
        static_assert(std::is_same_v<DFieldR, double> || ddc::is_borrowed_chunk_v<DFieldR>);
        if constexpr (ddc::is_borrowed_chunk_v<DFieldR>) {
            static_assert(DFieldR::rank() == 1);
        }
        return get_chunk_size(nustar);
    }

    // Get the size of the poloidal domain from B_norm by excluding the size of nustar
    template <class DFieldR, class DFieldRTheta>
    static std::size_t get_theta_domain_size(DFieldR nustar, DFieldRTheta B_norm)
    {
        std::size_t const r_size = get_chunk_size(nustar);
        std::size_t const r_theta_size = get_chunk_size(B_norm);
        return r_theta_size / r_size;
    }

    // Return a Kokkos view used in koliop from a DFieldR (used for the radial profile input)
    template <class DFieldR>
    koliop_interface::MDL<double*> radialprof_to_koliop(
            DFieldR radial_field,
            std::string const& string_radial_field)
    {
        koliop_interface::MDL<double*> kop_radial_field;

        std::size_t const n_r = get_r_domain_size(radial_field);
        if constexpr (std::is_same_v<DFieldR, double>) {
            kop_radial_field = koliop_interface::MDL<double*>(
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, string_radial_field),
                    n_r);
            Kokkos::deep_copy(kop_radial_field, radial_field);
        } else {
            kop_radial_field = radial_field.allocation_kokkos_view();
        }
        return kop_radial_field;
    }

    // Return a Kokkos view from B_norm of type DFieldRtheta. Enable to treat the case when n_r or/and ntheta equal to 1.
    template <class DFieldRTheta>
    koliop_interface::MDL<double**> rthetaprof_to_koliop(
            DFieldRTheta B_norm,
            std::size_t const n_r,
            std::size_t const n_theta)
    {
        koliop_interface::MDL<double**> kop_B_norm;

        if constexpr (std::is_same_v<DFieldRTheta, double>) {
            // If DFieldRTheta is a double create an array and copy the double to GPU
            kop_B_norm = koliop_interface::MDL<double**>(
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, "B_norm"),
                    n_r,
                    n_theta);
            Kokkos::deep_copy(kop_B_norm, B_norm);
        } else {
            // If DFieldRTheta is a chunk
            static_assert(
                    std::is_same_v<
                            typename DFieldRTheta::memory_space,
                            typename Kokkos::DefaultExecutionSpace::memory_space>,
                    "The data should be provided on GPU");
            auto B_norm_kokkos_view = B_norm.allocation_kokkos_view();
            if constexpr (DFieldRTheta::rank() == 1) {
                // If a 1D chunk then a reference to the data can be packed directly into a 2D chunk
                kop_B_norm = koliop_interface::MDL<double**>(B_norm.data_handle(), n_r, n_theta);
            } else if constexpr (
                    std::is_same_v<
                            typename decltype(B_norm_kokkos_view)::traits::array_layout,
                            typename koliop_interface::MDL<const double**>::traits::array_layout>) {
                // If a 2D chunk with the right layout we can use the kokkos view directly
                kop_B_norm = B_norm_kokkos_view;
            } else {
                // If a 2D chunk with the wrong layout then a copy is required
                static_assert(DFieldRTheta::rank() == 2);
                kop_B_norm = koliop_interface::MDL<double**>(
                        Kokkos::view_alloc(Kokkos::WithoutInitializing, "B_norm"),
                        n_r,
                        n_theta);
                // A copy is required to change the layout
                Kokkos::deep_copy(kop_B_norm, B_norm_kokkos_view);
            }
        }
        return kop_B_norm;
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
     * @param[in] nustar
     *      radial profile of collisionality
     * @param[in] rg
     *      radial profile of the grid. If size(nustar)==1, forced to 1. because already included in nustar definition 
     *      [TODO] See if this quantity cannot be included in nustar definition
     * @param[in] safety_factor
     *      radial profile of safety factor. If size(nustar)==1, forced to 1. because already included in nustar definition. 
     * @param[in] B_norm
     *      magnetic field norm in (r,theta)
     */
    template <class DFieldR, class DFieldRTheta>
    CollisionSpVparMu(
            FDistribDomain fdistrib_domain,
            DViewMu coeff_intdmu,
            DViewVpar coeff_intdvpar,
            DFieldR nustar,
            DFieldR rg, 
            DFieldR safety_factor,
            DFieldRTheta B_norm)
    : m_operator_handle {}
    , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
    , m_hat_As {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_As"), ddc::select<Species>(fdistrib_domain).size()}   // m_hat_As{"m_hat_As",ddc::select<Species>(fdistrib_domain)}
    , m_hat_Zs {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_hat_Zs"), ddc::select<Species>(fdistrib_domain).size()}
    , m_mask_buffer_r {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_buffer_r"), get_r_domain_size(nustar)}
    , m_mask_LIM {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mask_LIM"), get_r_domain_size(nustar), get_theta_domain_size(nustar, B_norm)}
    , m_Bstar_s {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_Bstar_s"),
              ddc::select<GridVpar>(fdistrib_domain).size(),
              get_r_domain_size(nustar),
              get_theta_domain_size(nustar, B_norm),
              ddc::select<Species>(fdistrib_domain).size()}
    , m_mug {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_mug"),
              ddc::select<GridMu>(fdistrib_domain).size()}
    , m_vparg {
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_vparg"),
              ddc::select<GridVpar>(fdistrib_domain).size()}
    {
        // Check that the distribution function is correctly ordered
        if constexpr (ddc::is_chunk_v<DFieldR>) {
            using GridR = domain_element_t<0, DFieldR>;
            if constexpr (ddc::is_chunk_v<DFieldRTheta>) {
                static_assert(std::is_same_v<GridR, domain_element_t<1, DFieldRTheta>>);
                using GridTheta = domain_element_t<0, DFieldRTheta>;
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
        } else if constexpr (ddc::is_chunk_v<DFieldRTheta>) {
            using GridTheta = domain_element_t<0, DFieldRTheta>;
            static_assert(
                    ddc::type_seq_rank_v<
                            GridTheta,
                            fdistrib_domain_tags> == (FDistribDomain::rank() - 3),
                    "Theta should appear third to last in the distribution function domain");
        }

        koliop_interface::DoCombMatComputation(m_comb_mat);

        IDomainSp idxrange_sp = ddc::select<Species>(fdistrib_domain);
        // --> Initialize the mass species
        ddc::ChunkSpan hat_As_host = ddc::discrete_space<IDimSp>().masses()[idxrange_sp];
        Kokkos::deep_copy(m_hat_As, hat_As_host.allocation_kokkos_view());
        // --> Initialize the charge species
        // TODO: Must be simplified as soon as the charge will be considered as a float and no more as an integer
        ddc::ChunkSpan hat_Zs_host = ddc::discrete_space<IDimSp>().charges()[idxrange_sp];
        auto hat_Zs = ddc::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                hat_Zs_host.span_cview());
        device_t<ddc::ChunkSpan<double, IDomainSp>> hat_Zs_proxy(m_hat_Zs.data(), idxrange_sp);
        ddc::parallel_deepcopy(Kokkos::DefaultExecutionSpace(), hat_Zs_proxy, hat_Zs);

        // --> Initialize the other quantities needed in koliop
        // TODO: Put Bstar_s as an input variable of the constructor (something more specific than what is done for B_norm must be done)
        Kokkos::deep_copy(m_Bstar_s, 1.0);

        // --> Initialization of the masks that have no sense here to 0.
        Kokkos::deep_copy(m_mask_buffer_r, 0.0); // Masked if >= 0.99
        Kokkos::deep_copy(m_mask_LIM, 0.0); // Masked if >= 0.99

        // --> Initialization of vpar and mu grids
        koliop_interface::DumpCoordinates(m_mug, ddc::select<GridMu>(fdistrib_domain));
        koliop_interface::DumpCoordinates(m_vparg, ddc::select<GridVpar>(fdistrib_domain));

        // NOTE: We need to fence because DumpCoordinates is asynchronous on GPUs.
        Kokkos::fence();

        std::size_t const n_mu = ddc::select<GridMu>(fdistrib_domain).size();
        std::size_t const n_vpar = ddc::select<GridVpar>(fdistrib_domain).size();
        std::size_t const n_r = get_r_domain_size(nustar);
        std::size_t const n_theta = get_theta_domain_size(nustar, B_norm);
        std::size_t const n_sp = ddc::select<Species>(fdistrib_domain).size();
        std::size_t const n_batch = fdistrib_domain.size() / (n_mu * n_vpar * n_r * n_theta * n_sp);

        m_nustar = radialprof_to_koliop(nustar, "nustar");
        m_rg = radialprof_to_koliop(rg, "rg");
        m_safety_factor = radialprof_to_koliop(safety_factor, "safety_factor");
        m_B_norm = rthetaprof_to_koliop(B_norm, n_r, n_theta);

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
                m_nustar.data(),
                m_comb_mat.data(),
                m_hat_As.data(), // m_hat_As.data_handle()   --> pointer
                m_hat_Zs.data(),
                m_rg.data(),
                m_safety_factor.data(),
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
    /// Radial profile of nustar
    koliop_interface::MDL<double*> m_nustar;
    /// Mesh points in the radial direction
    // [TODO]: See if we need m_rg ?
    koliop_interface::MDL<double*> m_rg;
    /// Radial safety factor profile
    koliop_interface::MDL<double*> m_safety_factor;
    /// Mask used to avoid to apply collision in certain region
    // [TODO]: This mask should maybe be deleted in C++ version
    koliop_interface::MDL<double*> m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    koliop_interface::MDL<double**> m_mask_LIM;
    /// B norm in (r,theta)
    // [TODO] Attention this must be 3D for generalization to 3D geometry--> transfer it in a 1D array ?
    koliop_interface::MDL<double**> m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // [TODO] Must be 5D for full 3D geometry
    koliop_interface::MDL<double****> m_Bstar_s;
    /// grid in mu direction
    koliop_interface::MDL<double*> m_mug;
    /// grid in vpar direction
    koliop_interface::MDL<double*> m_vparg;
};
