#pragma once

// Misc
#include "ddc/detail/tagged_vector.hpp"
#include "ddc/detail/type_seq.hpp"

// Containers
#include "ddc/aligned_allocator.hpp"
#include "ddc/chunk.hpp"
#include "ddc/chunk_span.hpp"

// Discretizations
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/discrete_vector.hpp"
#include "ddc/discretization.hpp"
#include "ddc/non_uniform_discretization.hpp"
#include "ddc/rectilinear_domain.hpp"
#include "ddc/uniform_discretization.hpp"
#include "ddc/uniform_domain.hpp"

// Algorithms
#include "ddc/deepcopy.hpp"
#include "ddc/for_each.hpp"
#include "ddc/reducer.hpp"
#include "ddc/transform_reduce.hpp"

// PDI wrapper
#if defined(DDC_BUILD_PDI_WRAPPER)
#include "ddc/pdi.hpp"
#endif
