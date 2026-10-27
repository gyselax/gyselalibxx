#!/usr/bin/env python3

import argparse
import sys

import h5py
import numpy as np


RELATIVE_ERROR_NORM2_TOLERANCE = 1e-10
MAX_ABSOLUTE_ERROR_TOLERANCE = 1e-10
MAX_RELATIVE_ERROR_TOLERANCE = 1e-10


def compare_hdf5(reference_file, approximation_file, dataset_name):
    with (
        h5py.File(reference_file, "r") as f_ref,
        h5py.File(approximation_file, "r") as f_approx,
    ):
        if dataset_name not in f_ref:
            raise KeyError(f"Dataset {dataset_name!r} not found in {reference_file}")
        if dataset_name not in f_approx:
            raise KeyError(
                f"Dataset {dataset_name!r} not found in {approximation_file}"
            )

        reference = np.asarray(f_ref[dataset_name][...], dtype=np.float64)
        approximation = np.asarray(f_approx[dataset_name][...], dtype=np.float64)

    if reference.shape != approximation.shape:
        raise ValueError(
            f"Shape mismatch: reference={reference.shape}, "
            f"approximation={approximation.shape}"
        )

    if reference.size == 0:
        raise ValueError("Cannot compare empty datasets")

    has_degenerated_reference = not np.all(np.isfinite(reference))
    has_degenerated_approximation = not np.all(np.isfinite(approximation))

    valid = np.isfinite(reference) & np.isfinite(approximation)

    if not np.any(valid):
        raise ValueError(
            "No entries are finite in both datasets; metrics are undefined"
        )

    r = reference[valid]
    a = approximation[valid]
    difference = a - r
    abs_error = np.abs(difference)
    abs_reference = np.abs(r)

    # L2 norms.
    reference_norm2 = np.linalg.norm(r)
    approximation_norm2 = np.linalg.norm(a)
    difference_norm2 = np.linalg.norm(difference)

    if reference_norm2 == 0.0:
        relative_error_norm2 = 0.0 if difference_norm2 == 0.0 else np.inf
    else:
        relative_error_norm2 = difference_norm2 / reference_norm2

    max_absolute_error = float(np.max(abs_error))

    relative_error = np.zeros_like(abs_error)
    nonzero_reference = abs_reference != 0.0

    np.divide(
        abs_error,
        abs_reference,
        out=relative_error,
        where=nonzero_reference,
    )
    relative_error[(~nonzero_reference) & (abs_error != 0.0)] = np.inf

    max_abs_flat = int(np.argmax(abs_error))
    max_rel_flat = int(np.argmax(relative_error))

    max_relative_error = float(relative_error[max_rel_flat])

    valid_indices = np.argwhere(valid)
    max_abs_index = tuple(valid_indices[max_abs_flat])
    max_rel_index = tuple(valid_indices[max_rel_flat])

    max_relative_error_at_max_absolute_error_location = float(
        relative_error[max_abs_flat]
    )
    max_absolute_error_at_max_relative_error_location = float(abs_error[max_rel_flat])

    warning = (
        has_degenerated_reference
        or has_degenerated_approximation
        or relative_error_norm2 > RELATIVE_ERROR_NORM2_TOLERANCE
        or (
            max_absolute_error > MAX_ABSOLUTE_ERROR_TOLERANCE
            and max_relative_error > MAX_RELATIVE_ERROR_TOLERANCE
        )
    )

    metrics = {
        "reference_norm2": reference_norm2,
        "approximation_norm2": approximation_norm2,
        "difference_norm2": difference_norm2,
        "relative_error_norm2": relative_error_norm2,
        "max_relative_error_at_max_absolute_error_location": max_relative_error_at_max_absolute_error_location,
        "max_relative_error_index": max_rel_index,
        "max_relative_error": max_relative_error,
        "max_absolute_error_at_max_relative_error_location": max_absolute_error_at_max_relative_error_location,
        "max_absolute_error_index": max_abs_index,
        "max_absolute_error": max_absolute_error,
        "has_degenerated_reference": has_degenerated_reference,
        "has_degenerated_approximation": has_degenerated_approximation,
        "warning": warning,
    }

    if warning:
        print("[WARNING]: Non conformity past tolerances detected!")

    name = dataset_name
    print(f"For dataset '{dataset_name}' with shape {reference.shape}:")
    print(
        f"NORM2 (r, a, a-r)           : {name} {reference_norm2:.17e} {approximation_norm2:.17e} {difference_norm2:.17e}"
    )
    print(f"NORM2(a-r) / NORM2(r)       : {name} {relative_error_norm2:.17e}")
    print(
        f"MAXVAL(ABS(a-r))            : {name} {max_absolute_error:.17e} ABS(a-r)/ABS(r) at argmax ABS(a-r): {max_relative_error_at_max_absolute_error_location:.17e}"
    )
    print(
        f"MAXVAL(ABS(a-r)/ABS(r))     : {name} {max_relative_error:.17e} ABS(a-r) at argmax relative error: {max_absolute_error_at_max_relative_error_location:.17e}"
    )
    print(f"Degenerated reference       : {name} {has_degenerated_reference}")
    print(f"Degenerated approximation   : {name} {has_degenerated_approximation}")
    print(f"Max absolute error index    : {max_abs_index}")
    print(f"Max relative error index    : {max_rel_index}")

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="Compare one dataset in two HDF5 files."
    )
    parser.add_argument("reference", help="Reference HDF5 file")
    parser.add_argument("approximation", help="Approximation HDF5 file")
    parser.add_argument("dataset", help="Dataset path inside both files")
    args = parser.parse_args()

    try:
        metrics = compare_hdf5(args.reference, args.approximation, args.dataset)
    except (OSError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 1 if metrics["warning"] else 0


if __name__ == "__main__":
    sys.exit(main())
