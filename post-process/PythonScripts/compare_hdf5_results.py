# SPDX-License-Identifier: MIT
""" Compare HDF5 results between two files.
"""
from argparse import ArgumentParser
import h5py as h5

import numpy as np

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Compare HDF5 results between two files')
    parser.add_argument('file1',
                        type=h5.File,
                        help='File name of the first HDF5 file')
    parser.add_argument('file2',
                        type=h5.File,
                        help='File name of the second HDF5 file')
    parser.add_argument('obj',
                        type=str,
                        help='Name of an HDF5 object, in absolute path')
    parser.add_argument('-R', '--relative',
            nargs='?',
            type=float,
            metavar='TOL',
            help='The tolerance for the relative difference')
    parser.add_argument('-A', '--absolute',
            nargs='?',
            type=float,
            metavar='TOL',
            help='The tolerance for the absolute difference')

    args = parser.parse_args()

    vals1 = np.array(args.file1[args.obj])
    vals2 = np.array(args.file2[args.obj])

    abs_error = np.abs(vals1 - vals2)
    rel_error = abs_error / np.abs(vals2)

    print("Maximum absolute error found : ", abs_error.max())
    print("Maximum relative error found : ", rel_error.max())

    args.file1.close()
    args.file2.close()

    if args.relative:
        assert args.absolute
        assert np.allclose(vals1, vals2, rtol=args.relative, atol=args.absolute)
    else:
        assert np.array_equal(vals1, vals2)
