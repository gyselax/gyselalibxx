// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

// Misc
#include "ddc/detail/macros.hpp"
#include "ddc/detail/tagged_vector.hpp"
#include "ddc/detail/type_seq.hpp"
#include "ddc/real_type.hpp"
#include "ddc/scope_guard.hpp"

// Containers
#include "ddc/aligned_allocator.hpp"
#include "ddc/chunk.hpp"
#include "ddc/chunk_span.hpp"
#include "ddc/kokkos_allocator.hpp"

// Discretizations
#include "ddc/coordinate_md.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_element.hpp"
#include "ddc/discrete_space.hpp"
#include "ddc/discrete_vector.hpp"
#include "ddc/non_uniform_point_sampling.hpp"
#include "ddc/periodic_sampling.hpp"
#include "ddc/rectilinear_domain.hpp"
#include "ddc/uniform_point_sampling.hpp"

// Algorithms
#include "ddc/for_each.hpp"
#include "ddc/mirror.hpp"
#include "ddc/parallel_deepcopy.hpp"
#include "ddc/parallel_fill.hpp"
#include "ddc/parallel_for_each.hpp"
#include "ddc/parallel_transform_reduce.hpp"
#include "ddc/reducer.hpp"
#include "ddc/transform_reduce.hpp"

// PDI wrapper
#if defined(DDC_BUILD_PDI_WRAPPER)
#include "ddc/pdi.hpp"
#endif

#if ginkgo_AVAIL
#include "misc/ginkgo_executors.hpp"
#endif
