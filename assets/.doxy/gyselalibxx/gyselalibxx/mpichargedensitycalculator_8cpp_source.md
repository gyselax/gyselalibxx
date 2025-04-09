

# File mpichargedensitycalculator.cpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**mpichargedensitycalculator.cpp**](mpichargedensitycalculator_8cpp.md)

[Go to the documentation of this file](mpichargedensitycalculator_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mpichargedensitycalculator.hpp"
#include "mpitools.hpp"
#include "species_info.hpp"

MpiChargeDensityCalculator::MpiChargeDensityCalculator(
        MPI_Comm comm,
        IChargeDensityCalculator const& local_charge_density_calculator)
    : m_local_charge_density_calculator(local_charge_density_calculator)
    , m_comm(comm)
{
}

void MpiChargeDensityCalculator::operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeDensityCalculator");

    DFieldMemXY rho_local_alloc(get_idx_range(rho));
    DFieldXY rho_local = get_field(rho_local_alloc);

    m_local_charge_density_calculator(rho_local, allfdistribu);

    Kokkos::DefaultExecutionSpace().fence("Fence local ChargeDensityCalculator");

    MPI_Allreduce(
            rho_local.data_handle(),
            rho.data_handle(),
            rho.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}
```


