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
        Chunk<cell, DiscreteDomain<DDimX, DDimY>>& cells)
{
    for_each(domain, [&](DiscreteCoordinate<DDimX, DDimY> const ixy) {
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
    DiscreteDomain<DDimX, DDimY> const domain_xy(
            DiscreteVector<DDimX, DDimY>(length, height));

    DiscreteDomain<DDimX, DDimY> const inner_domain_xy(
            DiscreteCoordinate<DDimX, DDimY>(1, 1),
            DiscreteVector<DDimX, DDimY>(length - 2, height - 2));

    Chunk<cell, DiscreteDomain<DDimX, DDimY>> cells_in(domain_xy);
    Chunk<cell, DiscreteDomain<DDimX, DDimY>> cells_out(domain_xy);

    // Initialize the whole domain
    blinker_init(domain_xy, cells_in);
    blinker_init(domain_xy, cells_out);

    std::size_t iter = 0;
    for (; iter < nt; ++iter) {
        print_2DChunk(std::cout, cells_in.span_cview()) << "\n";
        for_each(
                inner_domain_xy,
                [&](DiscreteCoordinate<DDimX, DDimY> const ixy) {
                    DiscreteCoordinate<DDimX> const ix
                            = select<DDimX>(ixy);
                    DiscreteCoordinate<DDimY> const iy
                            = select<DDimY>(ixy);
                    int alive_neighbors = 0;
                    // Iterate on neighbors and increase the count of alive neighbors when necessary
                    for (int i = -1; i < 2; ++i) {
                        for (int j = -1; j < 2; j++) {
                            if (cells_in(ix + i, iy + j)) {
                                std::cout << "ix : " << ix
                                          << " | i : " << i
                                          << " | ix + i : " << ix + i
                                          << "\n";
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
    print_2DChunk(std::cout, cells_in.span_cview()) << "\n";

    return 0;
}
