#pragma once

#include <ddc/ddc.hpp>

#include <assert.hpp>
#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "koliop_interface.hpp"

/**
 * @brief A namespace to collect classes which are necessary to create Chunks with the
 * correct number of dimensions to be compatible with Koliop.
 */
namespace collisions_dimensions {
/**
 * Create dimensions to act as the radial/poloidal/toroidal tag for Koliop even if these dims don't
 * exist in the simulation
 */

/// Class from which fake dimensions inherit. These fake dimensions are inserted in order to match Koliop interface
struct InternalSpoofGrid
{
};

/// Fake radial dimension to be used if there is no radial dimension in the simulation
struct InternalSpoofGridR : InternalSpoofGrid
{
};

/// Fake poloidal dimension to be used if there is no poloidal dimension in the simulation
struct InternalSpoofGridTheta : InternalSpoofGrid
{
};

/// Check if a dimension is spoofed but is not present in the actual simulation
template <class Grid>
inline constexpr bool is_spoofed_dim_v = std::is_base_of_v<InternalSpoofGrid, Grid>;

/**
 * Class to get the type of the radial dimension from a field containing a radial profile.
 * @tparam Field The type of the field containing the radial profile.
 */
template <class Field>
struct ExtractRDim
{
    static_assert(!std::is_same_v<Field, Field>, "Unrecognised radial profile type");
};

/**
 * Class to get the type of the poloidal dimension from a field containing a profile on the poloidal plane.
 * @tparam Field The type of the field containing the profile on the poloidal plane.
 * @tparam GridR The tag for the discrete radial dimension.
 */
template <class Field, class GridR>
struct ExtractThetaDim
{
    static_assert(!std::is_same_v<Field, Field>, "Unrecognised poloidal profile type");
};

/**
 * @brief Get the index range for a specific grid from a multi-D index range.
 * If the dimension is spoofed and does not appear in the multi-D index range then an index
 * range which only iterates over the index(0) is returned.
 *
 * @tparam Grid The tag for the specific grid.
 * @param idx_range The multi-D index range.
 * @returns The index range for the specific grid.
 */
template <class Grid, class FDistribDomain>
inline ddc::DiscreteDomain<Grid> get_1d_idx_range(FDistribDomain idx_range)
{
    if constexpr (is_spoofed_dim_v<Grid>) {
        return ddc::DiscreteDomain<
                Grid>(ddc::DiscreteElement<Grid> {0}, ddc::DiscreteVector<Grid> {1});
    } else {
        return ddc::select<Grid>(idx_range);
    }
}

/**
 * @brief Get the index range for specific grid dimensions from a multi-D index range.
 *
 * @tparam Grid The tags for the specific grid dimensions.
 * @param idx_range The multi-D index range.
 * @returns The index range for the specific grid.
 */
template <class... Grid, class FDistribDomain>
inline ddc::DiscreteDomain<Grid...> get_idx_range(FDistribDomain idx_range)
{
    return ddc::DiscreteDomain<Grid...>(get_1d_idx_range<Grid>(idx_range)...);
}

} // namespace collisions_dimensions

/**
 * @brief A class which computes the collision operator in (vpar,mu)
 */
template <
        class FDistribDomain,
        class GridVpar,
        class GridMu,
        class InputDFieldR,
        class InputDFieldRTheta>
class CollisionSpVparMu /* : public IRightHandSide */
{
private:
    using Species = IDimSp;
    using fdistrib_domain_tags = ddc::to_type_seq_t<FDistribDomain>;
    // Validate template types
    static_assert(ddc::is_discrete_domain_v<FDistribDomain>);
    static_assert(FDistribDomain::rank() >= 3 && FDistribDomain::rank() <= 6);
    static_assert((std::is_same_v<InputDFieldR, double>) || ddc::is_borrowed_chunk_v<InputDFieldR>);
    static_assert(
            (std::is_same_v<InputDFieldRTheta, double>)
            || (ddc::is_borrowed_chunk_v<InputDFieldRTheta>));
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

private:
    using GridR = typename collisions_dimensions::ExtractRDim<InputDFieldR>::type;
    using GridTheta =
            typename collisions_dimensions::ExtractThetaDim<InputDFieldRTheta, GridR>::type;

    using DDomR = ddc::DiscreteDomain<GridR>;
    using DDomRTheta = ddc::DiscreteDomain<GridR, GridTheta>;
    using DDomSpRThetaVpar = ddc::DiscreteDomain<GridVpar, GridR, GridTheta, Species>;

public:
    /// Type alias for the domain of the magnetic moment.
    using DDomMu = ddc::DiscreteDomain<GridMu>;
    /// Type alias for the domain of the velocity parallel to the magnetic field.
    using DDomVpar = ddc::DiscreteDomain<GridVpar>;

public:
    /// Type alias for a field on a grid of species
    using DFieldSp = device_t<ddc::Chunk<double, IDomainSp>>;
    /// Type alias for a field on a grid of radial values
    using DFieldR = device_t<ddc::Chunk<double, DDomR>>;
    /// Type alias for a field on a grid of magnetic moments
    using DFieldMu = device_t<ddc::Chunk<double, DDomMu>>;
    /// Type alias for a field on a grid of parallel velocities
    using DFieldVpar = device_t<ddc::Chunk<double, DDomVpar>>;
    /// Type alias for a field on a grid on a poloidal plane
    using DFieldRTheta = device_t<ddc::Chunk<double, DDomRTheta>>;
    /// Type alias for a field on a grid of species, poloidal plane and parallel velocities
    using DFieldSpRThetaVpar = device_t<ddc::Chunk<double, DDomSpRThetaVpar>>;
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
    template <class GridR, class SrcType>
    void deepcopy_1d(device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<GridR>>> dst, SrcType src)
    {
        if constexpr (collisions_dimensions::is_spoofed_dim_v<GridR>) {
            ddc::parallel_fill(dst, src);
        } else {
            ddc::parallel_deepcopy(dst, src);
        }
    }

    template <class GridR, class GridTheta, class SrcType>
    void deepcopy_2d(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<GridR, GridTheta>>> dst,
            SrcType src)
    {
        if constexpr ((collisions_dimensions::is_spoofed_dim_v<GridR>)&&(
                              collisions_dimensions::is_spoofed_dim_v<GridTheta>)) {
            ddc::parallel_fill(dst, src);
        } else if constexpr ((!collisions_dimensions::is_spoofed_dim_v<GridR>)&&(
                                     !collisions_dimensions::is_spoofed_dim_v<GridTheta>)) {
            ddc::parallel_deepcopy(dst, src);
        } else {
            using NonSpoofDim = std::
                    conditional_t<collisions_dimensions::is_spoofed_dim_v<GridR>, GridTheta, GridR>;
            ddc::parallel_for_each(dst.domain(), [&](ddc::DiscreteElement<GridR, GridTheta> idx) {
                dst(idx) = src(ddc::select<NonSpoofDim>(idx));
            });
        }
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
     * @param[in] collisions_interspecies
     *      boolean that is equal to true if inter-species collisions are taken into account
     * @param[in] rg
     *      radial profile of the grid. If size(nustar)==1, forced to 1. because already included in nustar definition 
     *      [TODO] See if this quantity cannot be included in nustar definition
     * @param[in] safety_factor
     *      radial profile of safety factor. If size(nustar)==1, forced to 1. because already included in nustar definition. 
     * @param[in] B_norm
     *      magnetic field norm in (r,theta)
     */
    CollisionSpVparMu(
            FDistribDomain fdistrib_domain,
            DViewMu coeff_intdmu,
            DViewVpar coeff_intdvpar,
            InputDFieldR nustar,
            std::int8_t const collisions_interspecies,
            InputDFieldR rg,
            InputDFieldR safety_factor,
            InputDFieldRTheta B_norm)
        : m_operator_handle {}
        , m_comb_mat {Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_comb_mat")}
        , m_hat_As {"m_hat_As", collisions_dimensions::get_idx_range<Species>(fdistrib_domain)}
        , m_hat_Zs {"m_hat_Zs", collisions_dimensions::get_idx_range<Species>(fdistrib_domain)}
        , m_nustar {"m_nustar", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_rg {"m_rg", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_safety_factor {"m_safety_factor", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_mask_buffer_r {"m_mask_buffer_r", collisions_dimensions::get_idx_range<GridR>(fdistrib_domain)}
        , m_mask_LIM {"m_mask_LIM", DDomRTheta {collisions_dimensions::get_idx_range<GridR, GridTheta>(fdistrib_domain)}}
        , m_B_norm {"m_B_norm", DDomRTheta {collisions_dimensions::get_idx_range<GridR, GridTheta>(fdistrib_domain)}}
        , m_Bstar_s {"m_Bstar_s", DDomSpRThetaVpar {collisions_dimensions::get_idx_range<Species, GridR, GridTheta, GridVpar>(fdistrib_domain)}}
        , m_mug {"m_mug", ddc::select<GridMu>(fdistrib_domain)}
        , m_vparg {"m_vparg", ddc::select<GridVpar>(fdistrib_domain)}
    {
        using namespace collisions_dimensions;
        // Check that the distribution function is correctly ordered
        if constexpr (!is_spoofed_dim_v<GridR>) {
            if constexpr (!is_spoofed_dim_v<GridTheta>) {
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
        } else if constexpr (!is_spoofed_dim_v<GridTheta>) {
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
        ddc::parallel_deepcopy(m_hat_As.span_view(), hat_As_host);
        // --> Initialize the charge species
        ddc::ChunkSpan hat_Zs_host = ddc::discrete_space<IDimSp>().charges()[idxrange_sp];
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
        std::size_t const n_r = get_idx_range<GridR>(fdistrib_domain).size();
        std::size_t const n_theta = get_idx_range<GridTheta>(fdistrib_domain).size();
        std::size_t const n_sp = ddc::select<Species>(fdistrib_domain).size();
        std::size_t const n_batch = fdistrib_domain.size() / (n_mu * n_vpar * n_r * n_theta * n_sp);

        deepcopy_1d<GridR>(m_nustar.span_view(), nustar);
        deepcopy_1d<GridR>(m_rg.span_view(), rg);
        deepcopy_1d<GridR>(m_safety_factor.span_view(), safety_factor);
        deepcopy_2d<GridR, GridTheta>(m_B_norm.span_view(), B_norm);

        m_operator_handle = koliop_interface::DoOperatorInitialization(
                n_mu,
                n_vpar,
                n_r,
                n_theta,
                n_batch,
                n_sp,
                collisions_interspecies,
                /* the_local_domain_r_offset */ 0 + n_r - 1,
                m_mug.data_handle(),
                m_vparg.data_handle(),
                coeff_intdmu.data_handle(),
                coeff_intdvpar.data_handle(),
                m_nustar.data_handle(),
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
    DFieldSp m_hat_As;
    /// Normalized charges for all species
    DFieldSp m_hat_Zs;
    /// Radial profile of nustar
    DFieldR m_nustar;
    /// Mesh points in the radial direction
    // [TODO]: See if we need m_rg ?
    DFieldR m_rg;
    /// Radial safety factor profile
    DFieldR m_safety_factor;
    /// Mask used to avoid to apply collision in certain region
    // [TODO]: This mask should maybe be deleted in C++ version
    DFieldR m_mask_buffer_r;
    /// Limiter mask in (r,theta)
    DFieldRTheta m_mask_LIM;
    /// B norm in (r,theta)
    // [TODO] Attention this must be 3D for generalization to 3D geometry--> transfer it in a 1D array ?
    DFieldRTheta m_B_norm;
    /// Bstar(species,r,theta,vpar)
    // [TODO] Must be 5D for full 3D geometry
    DFieldSpRThetaVpar m_Bstar_s;
    /// grid in mu direction
    DFieldMu m_mug;
    /// grid in vpar direction
    DFieldVpar m_vparg;
};

namespace collisions_dimensions {

/// If radial profile is stored in a double then the grid tag must be spoofed.
template <>
struct ExtractRDim<double>
{
    using type = InternalSpoofGridR;
};

/// If radial profile is stored in a 1D chunk then the grid tag is extracted.
template <class GridR, class Layout>
struct ExtractRDim<ddc::ChunkSpan<
        double,
        ddc::DiscreteDomain<GridR>,
        Layout,
        Kokkos::DefaultExecutionSpace::memory_space>>
{
    using type = GridR;
};

/// If the profile on the poloidal plane is stored in a double then the grid tag must be spoofed.
template <>
struct ExtractThetaDim<double, InternalSpoofGridR>
{
    using type = InternalSpoofGridTheta;
};

/// If the profile on the poloidal plane is stored in a chunk on radial values then the grid tag must be spoofed.
template <class GridR, class Layout>
struct ExtractThetaDim<
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<GridR>,
                Layout,
                Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = InternalSpoofGridTheta;
};

/// If the profile on the poloidal plane is stored in a chunk on poloidal values then the grid tag is extracted.
template <class GridR, class GridTheta, class Layout>
struct ExtractThetaDim<
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<GridTheta>,
                Layout,
                Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = GridTheta;
};

/// If the profile on the poloidal plane is stored in a 2D chunk then the grid tag is extracted.
template <class GridR, class GridTheta, class Layout>
struct ExtractThetaDim<
        ddc::ChunkSpan<
                double,
                ddc::DiscreteDomain<GridTheta, GridR>,
                Layout,
                Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = GridTheta;
};

} // namespace collisions_dimensions
