// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

using cell = bool;

// Name of the axis
struct DDimX;
struct DDimY;

static unsigned nt = 10;
static unsigned length = 5;
static unsigned height = 5;

void blinker_init(
        DiscreteDomain<DDimX, DDimY> const& domain,
        ChunkSpan<
                cell,
                DiscreteDomain<DDimX, DDimY>,
                std::experimental::layout_right,
                Kokkos::DefaultExecutionSpace::memory_space> cells)
{
    for_each(
            policies::parallel_device,
            domain,
            DDC_LAMBDA(DiscreteCoordinate<DDimX, DDimY> const ixy) {
                DiscreteCoordinate<DDimX> const ix = select<DDimX>(ixy);
                DiscreteCoordinate<DDimY> const iy = select<DDimY>(ixy);
                if (iy == DiscreteCoordinate<DDimY>(2)
                    && (ix >= DiscreteCoordinate<DDimX>(1)
                        && ix <= DiscreteCoordinate<DDimX>(3)))
                    cells(ixy) = true;
                else
                    cells(ixy) = false;
            });
}

template <class ElementType, class DDimX, class DDimY>
std::ostream& print_2DChunk(
        std::ostream& os,
        ChunkSpan<ElementType, DiscreteDomain<DDimX, DDimY>> chunk)
{
    for_each(
            select<DDimY>(chunk.domain()),
            [&](DiscreteCoordinate<DDimY> const iy) {
                for_each(
                        select<DDimX>(chunk.domain()),
                        [&](DiscreteCoordinate<DDimX> const ix) {
                            os << (chunk(ix, iy) ? "*" : ".");
                        });
                os << "\n";
            });
    return os;
}

int main()
{
    ScopeGuard scope;

    DiscreteDomain<DDimX, DDimY> const domain_xy(
            DiscreteCoordinate<DDimX, DDimY>(0, 0),
            DiscreteVector<DDimX, DDimY>(length, height));

    DiscreteDomain<DDimX, DDimY> const inner_domain_xy(
            DiscreteCoordinate<DDimX, DDimY>(1, 1),
            DiscreteVector<DDimX, DDimY>(length - 2, height - 2));

    Chunk cells_in_host_alloc(domain_xy, HostAllocator<cell>());
    Chunk cells_in_dev_alloc(domain_xy, DeviceAllocator<cell>());
    Chunk cells_out_dev_alloc(domain_xy, DeviceAllocator<cell>());

    ChunkSpan cells_in = cells_in_dev_alloc.span_view();
    ChunkSpan cells_out = cells_out_dev_alloc.span_view();

    // Initialize the whole domain
    blinker_init(domain_xy, cells_in);
    blinker_init(domain_xy, cells_out);

    std::size_t iter = 0;
    for (; iter < nt; ++iter) {
        deepcopy(cells_in_host_alloc, cells_in);
        print_2DChunk(std::cout, cells_in_host_alloc.span_cview())
                << "\n";
        for_each(
                policies::parallel_device,
                inner_domain_xy,
                DDC_LAMBDA(DiscreteCoordinate<DDimX, DDimY> const ixy) {
                    DiscreteCoordinate<DDimX> const ix
                            = select<DDimX>(ixy);
                    DiscreteCoordinate<DDimY> const iy
                            = select<DDimY>(ixy);
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
        deepcopy(cells_in, cells_out);
    }
    deepcopy(cells_in_host_alloc, cells_in);
    print_2DChunk(std::cout, cells_in_host_alloc.span_cview()) << "\n";

    return 0;
}
