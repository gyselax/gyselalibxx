#include <iomanip>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "sll/matrix_batch_tridiag.hpp"

#include "collisions_intra.hpp"
#include "collisions_utils.hpp"
#include "fluid_moments.hpp"

template <class TargetDim>
KOKKOS_FUNCTION Idx<TargetDim> CollisionsIntra::to_index(Idx<GridVx> const& index)
{
    static_assert(
            std::is_same_v<TargetDim, GhostedVx> || std::is_same_v<TargetDim, GhostedVxStaggered>);
    if constexpr (std::is_same_v<TargetDim, GhostedVx>) {
        return Idx<GhostedVx>(index.uid() + 1);
    } else {
        return Idx<GhostedVxStaggered>(index.uid() + 1);
    }
}

// The "Spoof" variables will be identical to the non-spoof versions. They are simply used
// to prevent the compiler from trying to compile code for the non-uniform case when splines
// are uniform.
template <class GridVxSpoof, class GhostedVxSpoof, class GhostedVxStaggeredSpoof>
std::enable_if_t<!ddc::is_uniform_point_sampling_v<GridVxSpoof>> CollisionsIntra::
        build_ghosted_staggered_vx_point_sampling(IdxRange<GridVxSpoof> const& idx_range)
{
    static_assert(
            std::is_same_v<GridVxSpoof, GridVx>,
            "The function is only designed to work with the GridVx dimension");
    static_assert(std::is_same_v<GhostedVxSpoof, GhostedVx>);
    static_assert(std::is_same_v<GhostedVxStaggeredSpoof, GhostedVxStaggered>);

    CoordVx const v0 = ddc::coordinate(idx_range.front());
    CoordVx const v1 = ddc::coordinate(idx_range.front() + 1);
    CoordVx const vN = ddc::coordinate(idx_range.back());
    CoordVx const vNm1 = ddc::coordinate(idx_range.back() - 1);
    int const ncells(idx_range.size() - 1);

    // ghosted points
    int const npoints(ncells + 3);
    std::vector<CoordVx> breaks(npoints);
    breaks[0] = v0 - (v1 - v0);
    breaks[npoints - 1] = vN + (vN - vNm1);
    ddc::for_each(idx_range, [&](IdxVx const iv) {
        breaks[to_index<GhostedVx>(iv).uid()] = ddc::coordinate(iv);
    });
    ddc::init_discrete_space<GhostedVxSpoof>(breaks);

    // ghosted staggered points
    int const npoints_stag(ncells + 2);
    std::vector<CoordVx> breaks_stag(npoints_stag);
    breaks_stag[0] = v0 - (v1 - v0) / 2.;
    breaks_stag[npoints_stag - 1] = vN + (vN - vNm1) / 2.;
    IdxRangeVx const gridv_less(idx_range.remove_last(IdxStepVx(1)));
    ddc::for_each(gridv_less, [&](IdxVx const iv) {
        breaks_stag[iv.uid() + 1] = CoordVx((ddc::coordinate(iv) + ddc::coordinate(iv + 1)) / 2.);
    });
    ddc::init_discrete_space<GhostedVxStaggeredSpoof>(breaks_stag);
}

// The "Spoof" variables will be identical to the non-spoof versions. They are simply used
// to prevent the compiler from trying to compile code for the uniform case when splines
// are non-uniform.
template <class GridVxSpoof, class GhostedVxSpoof, class GhostedVxStaggeredSpoof>
std::enable_if_t<ddc::is_uniform_point_sampling_v<GridVxSpoof>> CollisionsIntra::
        build_ghosted_staggered_vx_point_sampling(IdxRange<GridVxSpoof> const& idx_range)
{
    static_assert(
            std::is_same_v<GridVxSpoof, GridVx>,
            "The function is only designed to work with the GridVx dimension");
    static_assert(std::is_same_v<GhostedVxSpoof, GhostedVx>);
    static_assert(std::is_same_v<GhostedVxStaggeredSpoof, GhostedVxStaggered>);

    CoordVx const v0 = ddc::coordinate(idx_range.front());
    CoordVx const vN = ddc::coordinate(idx_range.back());
    int const ncells(idx_range.size() - 1);
    double const step(ddc::step<GridVxSpoof>());

    // ghosted points
    ddc::init_discrete_space<GhostedVx>(
            GhostedVxSpoof::init(v0 - step, vN + step, IdxStep<GhostedVx>(ncells + 3)));

    // ghosted staggered points
    ddc::init_discrete_space<GhostedVxStaggered>(
            GhostedVxStaggeredSpoof::
                    init(v0 - step / 2, vN + step / 2, IdxStep<GhostedVxStaggered>(ncells + 2)));
}

CollisionsIntra::CollisionsIntra(IdxRangeSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_fthresh(1.e-30)
    , m_nustar_profile_alloc(ddc::select<Species, GridX>(mesh))
    , m_gridvx_ghosted(Idx<GhostedVx>(0), IdxStep<GhostedVx>(ddc::select<GridVx>(mesh).size() + 2))
    , m_gridvx_ghosted_staggered(
              Idx<GhostedVxStaggered>(0),
              IdxStep<GhostedVxStaggered>(ddc::select<GridVx>(mesh).size() + 1))
    , m_mesh_ghosted(ddc::select<Species>(mesh), ddc::select<GridX>(mesh), m_gridvx_ghosted)
    , m_mesh_ghosted_staggered(
              ddc::select<Species>(mesh),
              ddc::select<GridX>(mesh),
              m_gridvx_ghosted_staggered)

{
    // validity checks
    if (ddc::select<Species>(mesh).size() != 2) {
        throw std::runtime_error("Intra species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    build_ghosted_staggered_vx_point_sampling(ddc::select<GridVx>(mesh));

    m_nustar_profile = get_field(m_nustar_profile_alloc);
    compute_nustar_profile(m_nustar_profile, m_nustar0);
    ddc::expose_to_pdi("collintra_nustar0", m_nustar0);
}

IdxRange<CollisionsIntra::GhostedVx> const& CollisionsIntra::get_gridvx_ghosted() const
{
    return m_gridvx_ghosted;
}

IdxRange<CollisionsIntra::GhostedVxStaggered> const& CollisionsIntra::get_gridvx_ghosted_staggered()
        const
{
    return m_gridvx_ghosted_staggered;
}

IdxRange<Species, GridX, CollisionsIntra::GhostedVx> const& CollisionsIntra::get_mesh_ghosted()
        const
{
    return m_mesh_ghosted;
}

void CollisionsIntra::compute_matrix_coeff(
        DFieldSpXVx AA,
        DFieldSpXVx BB,
        DFieldSpXVx CC,
        DField<IdxRangeSpXVx_ghosted> Dcoll,
        DField<IdxRangeSpXVx_ghosted_staggered> Dcoll_staggered,
        DField<IdxRangeSpXVx_ghosted> Nucoll,
        double deltat) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(AA),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);

                IdxVx_ghosted ivx_ghosted(to_index<GhostedVx>(ivx));
                IdxVx_ghosted_staggered ivx_ghosted_staggered(to_index<GhostedVxStaggered>(ivx));
                IdxVx_ghosted ivx_next_ghosted(ivx_ghosted + 1);
                IdxVx_ghosted ivx_prev_ghosted(ivx_ghosted - 1);
                IdxVx_ghosted_staggered ivx_prev_ghosted_staggered(ivx_ghosted_staggered - 1);

                double const dv_i
                        = ddc::coordinate(ivx_next_ghosted) - ddc::coordinate(ivx_ghosted);
                double const delta_i
                        = dv_i / (ddc::coordinate(ivx_ghosted) - ddc::coordinate(ivx_prev_ghosted));

                double const alpha_i = deltat / (dv_i * dv_i * (1. + delta_i));
                double const beta_i = deltat / (2. * dv_i * (1. + delta_i));

                double const coeffa
                        = alpha_i
                                  * (Dcoll_staggered(isp, ix, ivx_prev_ghosted_staggered) * delta_i
                                             * delta_i * delta_i
                                     - Dcoll(isp, ix, ivx_ghosted) * delta_i * delta_i
                                               * (delta_i - 1.))
                          + beta_i * Nucoll(isp, ix, ivx_prev_ghosted) * delta_i * delta_i;

                double const coeffb
                        = -alpha_i
                                  * (-Dcoll_staggered(isp, ix, ivx_ghosted_staggered)
                                     - Dcoll_staggered(isp, ix, ivx_prev_ghosted_staggered)
                                               * delta_i * delta_i * delta_i
                                     + Dcoll(isp, ix, ivx_ghosted) * (delta_i - 1.)
                                               * (delta_i * delta_i - 1.))
                          + beta_i * Nucoll(isp, ix, ivx_ghosted) * (delta_i * delta_i - 1.);

                double const coeffc = alpha_i
                                              * (Dcoll_staggered(isp, ix, ivx_ghosted_staggered)
                                                 + Dcoll(isp, ix, ivx_ghosted) * (delta_i - 1.))
                                      - beta_i * Nucoll(isp, ix, ivx_next_ghosted);

                AA(ispxvx) = -coeffa;
                BB(ispxvx) = 1. + coeffb;
                CC(ispxvx) = -coeffc;
            });
}

void CollisionsIntra::fill_matrix_with_coeff(
        Matrix_Banded& matrix,
        host_t<DConstFieldVx> AA,
        host_t<DConstFieldVx> BB,
        host_t<DConstFieldVx> CC) const
{
    matrix.set_element(0, 0, BB(IdxVx(0)));
    matrix.set_element(0, 1, CC(IdxVx(0)));

    int const npoints(get_idx_range<GridVx>(AA).size());
    matrix.set_element(npoints - 1, npoints - 1, BB(IdxVx(npoints - 1)));
    matrix.set_element(npoints - 1, npoints - 2, AA(IdxVx(npoints - 1)));

    IdxRangeVx const gridvx_inner(
            get_idx_range<GridVx>(AA).remove_first(IdxStepVx(1)).remove_last(IdxStepVx(1)));
    ddc::for_each(gridvx_inner, [&](IdxVx const ivx) {
        matrix.set_element(ivx.uid(), ivx.uid() - 1, AA(ivx));
        matrix.set_element(ivx.uid(), ivx.uid(), BB(ivx));
        matrix.set_element(ivx.uid(), ivx.uid() + 1, CC(ivx));
    });
}

void CollisionsIntra::compute_rhs_vector(
        DFieldSpXVx RR,
        DConstFieldSpXVx AA,
        DConstFieldSpXVx BB,
        DConstFieldSpXVx CC,
        DConstFieldSpXVx allfdistribu,
        double fthresh) const
{
    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(AA));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(RR),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);

                IdxVx const ivx_next = ivx + 1;
                IdxVx const ivx_prev = ivx - 1;

                if (ivx == idx_range_vx.front()) {
                    RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       + (-CC(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_next)
                                       - 2. * AA(isp, ix, ivx) * fthresh;

                } else if (ivx == idx_range_vx.back()) {
                    RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       + (-AA(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_prev)
                                       - 2. * CC(isp, ix, ivx) * fthresh;

                } else {
                    RR(isp, ix, ivx) = -AA(isp, ix, ivx) * allfdistribu(isp, ix, ivx_prev)
                                       + (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       - CC(isp, ix, ivx) * allfdistribu(isp, ix, ivx_next);
                }
            });
}



DFieldSpXVx CollisionsIntra::operator()(DFieldSpXVx allfdistribu, double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsIntra");

    IdxRangeSpX grid_sp_x(get_idx_range<Species, GridX>(allfdistribu));
    // density and temperature
    DFieldMemSpX density_alloc(grid_sp_x);
    DFieldMemSpX fluid_velocity_alloc(grid_sp_x);
    DFieldMemSpX temperature_alloc(grid_sp_x);
    DFieldSpX density = get_field(density_alloc);
    DFieldSpX fluid_velocity = get_field(fluid_velocity_alloc);
    DFieldSpX temperature = get_field(temperature_alloc);

    DFieldMemVx const quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(
                    get_idx_range<GridVx>(allfdistribu)));
    DConstFieldVx quadrature_coeffs = get_const_field(quadrature_coeffs_alloc);

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    //Moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid_sp_x,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : idx_range_vx) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                fluid_velocity(ispx) = particle_flux / density(ispx);
                temperature(ispx)
                        = (momentum_flux - particle_flux * fluid_velocity(ispx)) / density(ispx);
            });

    // collision frequency
    DFieldMemSpX collfreq_alloc(grid_sp_x);
    DFieldSpX collfreq = get_field(collfreq_alloc);
    compute_collfreq(
            collfreq,
            get_const_field(m_nustar_profile),
            get_const_field(density),
            get_const_field(temperature));

    // diffusion coefficient
    DFieldMem<IdxRangeSpXVx_ghosted> Dcoll_alloc(m_mesh_ghosted);
    DField<IdxRangeSpXVx_ghosted> Dcoll = get_field(Dcoll_alloc);
    compute_Dcoll<GhostedVx>(
            Dcoll,
            get_const_field(collfreq),
            get_const_field(density),
            get_const_field(temperature));

    DFieldMem<IdxRangeSpXVx_ghosted> dvDcoll_alloc(m_mesh_ghosted);
    DField<IdxRangeSpXVx_ghosted> dvDcoll = get_field(dvDcoll_alloc);
    compute_dvDcoll<GhostedVx>(
            dvDcoll,
            get_const_field(collfreq),
            get_const_field(density),
            get_const_field(temperature));

    DFieldMem<IdxRangeSpXVx_ghosted_staggered> Dcoll_staggered_alloc(m_mesh_ghosted_staggered);
    DField<IdxRangeSpXVx_ghosted_staggered> Dcoll_staggered = get_field(Dcoll_staggered_alloc);
    compute_Dcoll<GhostedVxStaggered>(
            Dcoll_staggered,
            get_const_field(collfreq),
            get_const_field(density),
            get_const_field(temperature));

    // kernel maxwellian fluid moments
    DFieldMemSpX Vcoll_alloc(grid_sp_x);
    DFieldMemSpX Tcoll_alloc(grid_sp_x);
    DFieldSpX Vcoll = get_field(Vcoll_alloc);
    DFieldSpX Tcoll = get_field(Tcoll_alloc);
    compute_Vcoll_Tcoll<GhostedVx>(Vcoll, Tcoll, get_const_field(allfdistribu), Dcoll, dvDcoll);

    // convection coefficient Nucoll
    DFieldMem<IdxRangeSpXVx_ghosted> Nucoll_alloc(m_mesh_ghosted);
    DField<IdxRangeSpXVx_ghosted> Nucoll = get_field(Nucoll_alloc);
    compute_Nucoll<GhostedVx>(Nucoll, Dcoll, get_const_field(Vcoll), get_const_field(Tcoll));

    // matrix coefficients
    DFieldMemSpXVx AA_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx BB_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx CC_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx AA = get_field(AA_alloc);
    DFieldSpXVx BB = get_field(BB_alloc);
    DFieldSpXVx CC = get_field(CC_alloc);
    compute_matrix_coeff(AA, BB, CC, Dcoll, Dcoll_staggered, Nucoll, dt);


    // rhs vector coefficient
    DFieldMemSpXVx RR_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx RR = get_field(RR_alloc);
    compute_rhs_vector(
            RR,
            get_const_field(AA),
            get_const_field(BB),
            get_const_field(CC),
            get_const_field(allfdistribu),
            m_fthresh);

    int const batch_size = get_idx_range<Species, GridX>(allfdistribu).size();
    int const mat_size = get_idx_range<GridVx>(allfdistribu).size();
    /* Here we do not use allocation_kokkos_view() ddc function since we change the shape 
       from (Sp,X,Vx)-->(batch_dim,Vx)*/
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            AA_view(AA.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            BB_view(BB.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            CC_view(CC.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            RR_view(RR.data_handle(), batch_size, mat_size);


    MatrixBatchTridiag<Kokkos::DefaultExecutionSpace>
            matrix(batch_size, mat_size, AA_view, BB_view, CC_view);
    matrix.setup_solver();
    matrix.solve(RR_view);

    ddc::parallel_deepcopy(allfdistribu, RR);
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
