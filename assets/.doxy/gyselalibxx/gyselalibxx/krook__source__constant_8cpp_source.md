

# File krook\_source\_constant.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**krook\_source\_constant.cpp**](krook__source__constant_8cpp.md)

[Go to the documentation of this file](krook__source__constant_8cpp.md)


```C++
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "krook_source_constant.hpp"
#include "mask_tanh.hpp"
#include "maxwellianequilibrium.hpp"

KrookSourceConstant::KrookSourceConstant(
        IdxRangeX const& gridx,
        IdxRangeVx const& gridvx,
        RhsType const type,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const temperature)
    : m_type(type)
    , m_extent(extent)
    , m_stiffness(stiffness)
    , m_amplitude(amplitude)
    , m_density(density)
    , m_temperature(temperature)
    , m_mask(gridx)
    , m_ftarget(gridvx)
{
    // mask that defines the region where the operator is active
    host_t<DFieldMemX> mask_host(gridx);
    switch (m_type) {
    case RhsType::Source:
        // the mask equals one in the interval [x_left, x_right]
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Normal, false);
        break;
    case RhsType::Sink:
        // the mask equals zero in the centre of the plasma
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Inverted, false);
        break;
    }
    ddc::parallel_deepcopy(get_field(m_mask), mask_host);

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(get_field(m_ftarget), m_density, m_temperature, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(get_field(m_ftarget));

    switch (m_type) {
    case RhsType::Source:
        ddc::expose_to_pdi("krook_source_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_source_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_constant_density", m_density);
        ddc::expose_to_pdi("krook_source_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_constant_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_source_constant_mask", mask_host);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_constant_density", m_density);
        ddc::expose_to_pdi("krook_sink_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_constant_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_sink_constant_mask", mask_host);
        break;
    }
}

DFieldSpXVx KrookSourceConstant::operator()(DFieldSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("KrookSource");

    DConstFieldVx ftarget = get_const_field(m_ftarget);
    DConstFieldX mask = get_const_field(m_mask);
    double const amplitude = m_amplitude;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                allfdistribu(ispxvx)
                        = ftarget(ddc::select<GridVx>(ispxvx))
                          + (allfdistribu(ispxvx) - ftarget(ddc::select<GridVx>(ispxvx)))
                                    * Kokkos::exp(
                                            -amplitude * mask(ddc::select<GridX>(ispxvx)) * dt);
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
```


