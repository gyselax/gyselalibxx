# Copyright (C) The DDC development team, see COPYRIGHT.md file
#
# SPDX-License-Identifier: MIT

# FindLAPACKE.cmake
# ----------------
# Finds the LAPACKE library (C interface to LAPACK).
#
# LAPACKE may be bundled inside an existing LAPACK/BLAS provider (e.g.
# OpenBLAS) or shipped as an independent liblapacke. This module first
# delegates to CMake's own FindLAPACK to locate the LAPACK libraries, then
# checks whether LAPACKE symbols are already present in those libraries
# before falling back to searching for a standalone liblapacke.
#
# The header lapacke.h is always resolved from the shared object that actually
# provides the LAPACKE symbols, so that the include path and the library stay
# in sync.
#
# Imported target
# ^^^^^^^^^^^^^^^
#   LAPACKE::LAPACKE
#     Compile-and-link target. Carries the correct include directory and links
#     to whichever library (bundled or standalone) provides LAPACKE.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#   LAPACKE_FOUND          - TRUE when both the library and the header are found
#   LAPACKE_LIBRARIES      - Full list of libraries needed to use LAPACKE
#                            (may re-use LAPACK_LIBRARIES when bundled)
#   LAPACKE_INCLUDE_DIRS   - Include directory containing lapacke.h
#
# Hints
# ^^^^^
#   LAPACKE_ROOT / ENV{LAPACKE_ROOT}
#     Prefix searched first for the standalone library and header.
#   LAPACK_ROOT / ENV{LAPACK_ROOT}  (also honoured by FindLAPACK)
#     Prefix for the underlying LAPACK provider.
#
# Notes
# ^^^^^
# * BLA_VENDOR (see FindLAPACK docs) can be set before find_package(LAPACKE)
#   to steer which BLAS/LAPACK flavour is considered, e.g.
#     set(BLA_VENDOR OpenBLAS)
#     find_package(LAPACKE REQUIRED)
# * The check_function_exists() probes require try_compile, so the C language
#   must be enabled in the calling project.

# Early exit: if a previous call already succeeded, do not redo all the work.
# LAPACKE_FOUND is cached by find_package_handle_standard_args on success.
if(LAPACKE_FOUND)
    return()
endif()

if(CMAKE_VERSION VERSION_LESS "3.25")
    message(FATAL_ERROR "CMake >= 3.25 required")
endif()

include(CheckFunctionExists)
include(FindPackageHandleStandardArgs)

# ---------------------------------------------------------------------------
# Step 1 – Locate LAPACK (and its transitive BLAS dependency).
# ---------------------------------------------------------------------------
find_package(LAPACK QUIET)

# ---------------------------------------------------------------------------
# Step 2 – Check whether LAPACKE is already bundled in the LAPACK libraries.
#
# OpenBLAS and MKL both ship LAPACKE symbols inside their main library, so
# we avoid adding a redundant -llapacke in those cases.
# ---------------------------------------------------------------------------
set(_lapacke_bundled FALSE)

if(LAPACK_FOUND)
    # Save/restore CMAKE_REQUIRED_LIBRARIES to avoid polluting the caller's state.
    # Only LAPACK_LIBRARIES is needed: a correctly built shared library encodes
    # its BLAS dependency as a DT_NEEDED entry, so the dynamic linker resolves
    # it without an explicit -lblas flag. Passing BLAS_LIBRARIES here would mask
    # a liblapack with missing DT_NEEDED entries.
    set(_saved_req_libs "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})

    # LAPACKE_dgetrf (LU factorisation) is a fundamental routine present in
    # every complete LAPACKE implementation and a reliable probe symbol.
    check_function_exists(LAPACKE_dgetrf _lapacke_bundled_check)

    set(CMAKE_REQUIRED_LIBRARIES "${_saved_req_libs}")
    unset(_saved_req_libs)

    if($CACHE{_lapacke_bundled_check})
        set(_lapacke_bundled TRUE)
        message(STATUS "FindLAPACKE: LAPACKE symbols found inside LAPACK libraries (bundled)")
    endif()

endif()

# ---------------------------------------------------------------------------
# Step 3 – If not bundled, search for a standalone liblapacke.
# ---------------------------------------------------------------------------
if(NOT _lapacke_bundled)
    find_library(LAPACKE_LIBRARY
        NAMES lapacke
        HINTS
            "${LAPACKE_ROOT}"
            ENV LAPACKE_ROOT
            "${LAPACK_ROOT}"
            ENV LAPACK_ROOT
        PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu lib/aarch64-linux-gnu
    )

    if(LAPACKE_LIBRARY)
        # Verify the symbols are really there (guards against empty stub libs).
        # Only liblapacke itself is needed: a correctly built shared library
        # encodes its LAPACK/BLAS dependencies as DT_NEEDED entries, so the
        # dynamic linker resolves them without explicit -llapack/-lblas flags.
        # Passing them here would mask a liblapacke with missing DT_NEEDED.
        set(_saved_req_libs "${CMAKE_REQUIRED_LIBRARIES}")
        set(CMAKE_REQUIRED_LIBRARIES "${LAPACKE_LIBRARY}")

        check_function_exists(LAPACKE_dgetrf _lapacke_standalone_check)

        set(CMAKE_REQUIRED_LIBRARIES "${_saved_req_libs}")
        unset(_saved_req_libs)

        if($CACHE{_lapacke_standalone_check})
            message(STATUS "FindLAPACKE: found standalone LAPACKE library: ${LAPACKE_LIBRARY}")
        else()
            message(WARNING "FindLAPACKE: liblapacke found at ${LAPACKE_LIBRARY} "
                            "but LAPACKE_dgetrf is not resolvable — ignoring.")
            unset(LAPACKE_LIBRARY CACHE)
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# Step 4 – Locate lapacke.h from the library that provides the symbols.
# ---------------------------------------------------------------------------
# Prefer an explicit hint, then probe near the providing library.
find_path(LAPACKE_INCLUDE_DIR
    NAMES lapacke.h
    HINTS
        "${LAPACKE_ROOT}"
        ENV LAPACKE_ROOT
        "${LAPACK_ROOT}"
        ENV LAPACK_ROOT
    PATH_SUFFIXES include include/openblas include/lapacke include/lapack
)

# ---------------------------------------------------------------------------
# Step 5 – Assemble result variables.
# ---------------------------------------------------------------------------
if(_lapacke_bundled)
    # All required symbols live in LAPACK_LIBRARIES; no extra -l needed.
    set(LAPACKE_LIBRARIES ${LAPACK_LIBRARIES})
elseif(LAPACKE_LIBRARY)
    # Standalone: prepend liblapacke, keep LAPACK/BLAS for transitive deps.
    set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARIES})
else()
    set(LAPACKE_LIBRARIES "")
endif()

set(LAPACKE_INCLUDE_DIRS "${LAPACKE_INCLUDE_DIR}")

# ---------------------------------------------------------------------------
# Step 6 – Standard "found" handling.
# ---------------------------------------------------------------------------
find_package_handle_standard_args(LAPACKE
    REQUIRED_VARS
        LAPACKE_LIBRARIES
        LAPACKE_INCLUDE_DIRS
    REASON_FAILURE_MESSAGE
        "Could not find LAPACKE. Set LAPACKE_ROOT or BLA_VENDOR to guide the search."
)

# ---------------------------------------------------------------------------
# Step 7 – Create the imported target.
# ---------------------------------------------------------------------------
if(LAPACKE_FOUND AND NOT TARGET LAPACKE::LAPACKE)
    add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
    set_target_properties(LAPACKE::LAPACKE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES      "${LAPACKE_LIBRARIES}"
    )
endif()

unset(_lapacke_bundled)
