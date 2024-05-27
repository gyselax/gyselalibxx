// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

using cell = bool;

// Name of the axis
struct DDimX
{
};
struct DDimY
{
};

static unsigned nt = 10;
static unsigned length = 5;
static unsigned height = 5;

void blinker_init(
        ddc::DiscreteDomain<DDimX, DDimY> const& domain,
        ddc::ChunkSpan<
                cell,
                ddc::DiscreteDomain<DDimX, DDimY>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> cells)
{
    ddc::parallel_for_each(
            domain,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                ddc::DiscreteElement<DDimX> const ix
                        = ddc::select<DDimX>(ixy);
                ddc::DiscreteElement<DDimY> const iy
                        = ddc::select<DDimY>(ixy);
                if (iy == ddc::DiscreteElement<DDimY>(2)
                    && (ix >= ddc::DiscreteElement<DDimX>(1)
                        && ix <= ddc::DiscreteElement<DDimX>(3)))
                    cells(ixy) = true;
                else
                    cells(ixy) = false;
            });
}

template <class ElementType, class DDimX, class DDimY>
std::ostream& print_2DChunk(
        std::ostream& os,
        ddc::ChunkSpan<ElementType, ddc::DiscreteDomain<DDimX, DDimY>>
                chunk)
{
    ddc::for_each(
            ddc::select<DDimY>(chunk.domain()),
            [&](ddc::DiscreteElement<DDimY> const iy) {
                ddc::for_each(
                        ddc::select<DDimX>(chunk.domain()),
                        [&](ddc::DiscreteElement<DDimX> const ix) {
                            os << (chunk(ix, iy) ? "*" : ".");
                        });
                os << "\n";
            });
    return os;
}

int main()
{
    Kokkos::ScopeGuard const kokkos_scope;
    ddc::ScopeGuard const ddc_scope;

    ddc::DiscreteDomain<DDimX, DDimY> const domain_xy(
            ddc::DiscreteElement<DDimX, DDimY>(0, 0),
            ddc::DiscreteVector<DDimX, DDimY>(length, height));

    ddc::DiscreteDomain<DDimX, DDimY> const inner_domain_xy(
            ddc::DiscreteElement<DDimX, DDimY>(1, 1),
            ddc::DiscreteVector<DDimX, DDimY>(length - 2, height - 2));

    ddc::Chunk cells_in_host_alloc(
            "cells_in_host",
            domain_xy,
            ddc::HostAllocator<cell>());
    ddc::Chunk cells_in_dev_alloc(
            "cells_in_dev",
            domain_xy,
            ddc::DeviceAllocator<cell>());
    ddc::Chunk cells_out_dev_alloc(
            "cells_out_dev",
            domain_xy,
            ddc::DeviceAllocator<cell>());

    ddc::ChunkSpan cells_in = cells_in_dev_alloc.span_view();
    ddc::ChunkSpan cells_out = cells_out_dev_alloc.span_view();

    // Initialize the whole domain
    blinker_init(domain_xy, cells_in);
    blinker_init(domain_xy, cells_out);

    std::size_t iter = 0;
    for (; iter < nt; ++iter) {
        ddc::parallel_deepcopy(cells_in_host_alloc, cells_in);
        print_2DChunk(std::cout, cells_in_host_alloc.span_cview())
                << "\n";
        ddc::parallel_for_each(
                inner_domain_xy,
                KOKKOS_LAMBDA(
                        ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                    ddc::DiscreteElement<DDimX> const ix
                            = ddc::select<DDimX>(ixy);
                    ddc::DiscreteElement<DDimY> const iy
                            = ddc::select<DDimY>(ixy);
                    int alive_neighbors = 0;
                    // Iterate on neighbors and increase the count of alive neighbors when necessary
                    for (int i = -1; i < 2; ++i) {
                        for (int j = -1; j < 2; j++) {
                            if (cells_in(ix + i, iy + j)) {
                                alive_neighbors++;
                            }
                        }
                    }
                    // Update the future status of the current cell depending on its current status and
                    // its current number of alive neighbors
                    if (cells_in(ixy)) {
                        alive_neighbors--;
                        if (alive_neighbors < 2 || alive_neighbors > 3)
                            cells_out(ixy) = false;
                    } else {
                        if (alive_neighbors == 3)
                            cells_out(ixy) = true;
                    }
                });
        ddc::parallel_deepcopy(cells_in, cells_out);
    }
    ddc::parallel_deepcopy(cells_in_host_alloc, cells_in);
    print_2DChunk(std::cout, cells_in_host_alloc.span_cview()) << "\n";

    return 0;
}
