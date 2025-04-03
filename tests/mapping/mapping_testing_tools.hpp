// SPDX-License-Identifier: MIT
#pragma once

template <class GridR, class GridTheta>
auto get_example_coords(IdxStep<GridR> r_ncells, IdxStep<GridTheta> theta_ncells)
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;

    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxRTheta = Idx<GridR, GridTheta>;

    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;

    Coord<R> const r_min(0.0);
    Coord<R> const r_max(1.0);

    Coord<Theta> const theta_min(0.0);
    Coord<Theta> const theta_max(2.0 * M_PI);

    IdxR const r_zero(0);
    IdxR const r_start(1); // avoid singular point at r = 0.
    IdxTheta const theta_start(0);

    double const dr((r_max - r_min) / r_ncells);
    double const dtheta((theta_max - theta_min) / theta_ncells);

    IdxRangeR idx_range_r(r_start, r_ncells);
    IdxRangeTheta idx_range_theta(theta_start, theta_ncells);
    IdxRangeRTheta grid(idx_range_r, idx_range_theta);

    host_t<FieldMem<Coord<R, Theta>, IdxRangeRTheta>> coords(grid);
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        IdxR ir(irtheta);
        IdxTheta itheta(irtheta);
        Coord<R> r(r_min + dr * (ir - r_zero));
        Coord<Theta> theta(theta_min + dtheta * (itheta - theta_start));
        coords(irtheta) = Coord<R, Theta>(r, theta);
    });

    return coords;
}

template <class GridR, class GridTheta, class GridPhi>
auto get_example_coords(
        IdxStep<GridR> r_ncells,
        IdxStep<GridTheta> theta_ncells,
        IdxStep<GridPhi> phi_ncells)
{
    using R = typename GridR::continuous_dimension_type;
    using Theta = typename GridTheta::continuous_dimension_type;
    using Phi = typename GridPhi::continuous_dimension_type;

    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxPhi = Idx<GridPhi>;
    using IdxRThetaPhi = Idx<GridR, GridTheta, GridPhi>;

    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangePhi = IdxRange<GridPhi>;
    using IdxRangeRThetaPhi = IdxRange<GridR, GridTheta, GridPhi>;

    Coord<R> const r_min(0.0);
    Coord<R> const r_max(1.0);

    Coord<Theta> const theta_min(0.0);
    Coord<Theta> const theta_max(2.0 * M_PI);

    Coord<Phi> const phi_min(0.0);
    Coord<Phi> const phi_max(2.0 * M_PI);

    IdxR const r_zero(0);
    IdxR const r_start(1); // avoid singular point at r = 0.
    IdxTheta const theta_start(0);
    IdxPhi const phi_start(0);

    double const dr((r_max - r_min) / r_ncells);
    double const dtheta((theta_max - theta_min) / theta_ncells);
    double const dphi((phi_max - phi_min) / phi_ncells);

    IdxRangeR idx_range_r(r_start, r_ncells);
    IdxRangeTheta idx_range_theta(theta_start, theta_ncells);
    IdxRangePhi idx_range_phi(phi_start, phi_ncells);
    IdxRangeRThetaPhi grid(idx_range_r, idx_range_theta, idx_range_phi);

    host_t<FieldMem<Coord<R, Theta, Phi>, IdxRangeRThetaPhi>> coords(grid);
    ddc::for_each(grid, [&](IdxRThetaPhi const irthetaphi) {
        IdxR ir(irthetaphi);
        IdxTheta itheta(irthetaphi);
        IdxPhi iphi(irthetaphi);
        Coord<R> r(r_min + dr * (ir - r_zero));
        Coord<Theta> theta(theta_min + dtheta * (itheta - theta_start));
        Coord<Phi> phi(phi_min + dphi * (iphi - phi_start));
        coords(irthetaphi) = Coord<R, Theta, Phi>(r, theta, phi);
    });

    return coords;
}
