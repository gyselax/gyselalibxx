import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import h5py

description="""Script to provide break points for non-equidistant points.
This script allows the user to specify domains and the number of cells in those domains

The domains are specified with the edge-domains parameter.
The small domains are defined between consecutive edge-domains.
In the periodic case, the final domain wraps from the last edge point, round to the first.
So the domains are defined as follows:
    [e_1, e_2] , [e_2, e_3], ... , [e_{n-1}, e_n] , [e_n, e_1]

For each domain the number of cells must be specified with the n-cells parameter.
Therefore the number of cells must be equal to the number of edge domains.
"""

def get_non_periodic_points(n_cells, edge_domains):
    if len(edge_domains) < 2:
        print("ERROR : There must be at least one domain which must contain at least 1 cell.")
        sys.exit(1)
    if len(n_cells) != len(edge_domains) - 1:
        print("ERROR : There are", len(edge_domains) - 1, "domains, but number of cells were provided for", len(n_cells), "domains")
        sys.exit(1)

    n_domains = len(n_cells)

    points = [np.linspace(edge_domains[i], edge_domains[i+1], n_cells[i], endpoint = False) for i in range(n_domains)]
    return [p for dom in points for p in dom] + [edge_domains[-1]]

def get_periodic_points(xmin, xmax, n_cells, edge_domains):
    if edge_domains[0] == xmin:
        return get_non_periodic_points(n_cells, edge_domains)[:-1]
    if len(edge_domains) < 1:
        print("ERROR : There must be at least one domain which must contain at least 1 cell.")
        sys.exit(1)
    if len(n_cells) != len(edge_domains):
        print("ERROR : There are", len(edge_domains), "domains, but number of cells were provided for", len(n_cells), "domains")
        sys.exit(1)

    n_domains = len(n_cells)
    Lx = xmax - xmin

    edges = [xmin, *edge_domains, xmax]
    wrapped_dom_len = edge_domains[0]+Lx - edge_domains[-1]
    n_wrapped = n_cells[-1]
    n_split_1 = max(1, int(n_wrapped*(edge_domains[0]-xmin)/wrapped_dom_len))
    n_split_2 = n_wrapped - n_split_1
    points = [n_split_1, *n_cells[:-1], n_split_2]
    assert sum(points) == sum(n_cells)

    points = [np.linspace(edges[i], edges[i+1], points[i], endpoint = False) for i in range(n_domains+1)]
    return [p for dom in points for p in dom] + [xmax]

if __name__ == '__main__':
    parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('out_file', type=str, help="Name of the hdf5 output file")
    parser.add_argument('--name', type=str, help="The name under which the new grid should be saved", required=True)
    parser.add_argument('--edge-domains', nargs='+', type=float, help="Coordinates of the edge of the domains", required=True)
    parser.add_argument('--ncells', nargs='+', type=int, help="Number of cells in each domain", required=True)
    parser.add_argument('--xmin', type=float, help="Start coordinate of the full domain", required=False)
    parser.add_argument('--xmax', type=float, help="End coordinate of the full domain", required=False)
    parser.add_argument('--periodic', action='store_true', help="Indicate that the domain is periodic")
    args = parser.parse_args()

    edge_domains = args.edge_domains
    n_cells = args.ncells
    xmin = args.xmin
    xmax = args.xmax

    if args.periodic:
        if xmin is None:
            raise RuntimeError("xmin must be provided in the periodic case")
        if xmax is None:
            raise RuntimeError("xmax must be provided in the periodic case")
        points = get_periodic_points(xmin, xmax, n_cells, edge_domains)
    else:
        points = get_non_periodic_points(n_cells, edge_domains)

    print(points, len(points))

    h5 = h5py.File( args.out_file, mode='a' )
    h5.create_dataset( args.name, data=points )
    h5.close()

