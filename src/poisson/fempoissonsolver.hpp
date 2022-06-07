#pragma once

#if defined(ENABLE_PERIODIC_RDIMX)
#if ENABLE_PERIODIC_RDIMX

#include "femperiodicpoissonsolver.hpp"

using FemPoissonSolver = FemPeriodicPoissonSolver;

#else

#include "femnonperiodicpoissonsolver.hpp"

using FemPoissonSolver = FemNonPeriodicPoissonSolver;

#endif
#endif
