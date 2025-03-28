

# File fft\_poisson\_solver.hpp

[**File List**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**fft\_poisson\_solver.hpp**](fft__poisson__solver_8hpp.md)

[Go to the documentation of this file](fft__poisson__solver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "ipoisson_solver.hpp"
#include "vector_index_tools.hpp"

template <
        class IdxRangeLaplacian,
        class IdxRangeFull = IdxRangeLaplacian,
        class ExecSpace = Kokkos::DefaultExecutionSpace,
        class LayoutSpace = Kokkos::layout_right>
class FFTPoissonSolver;

template <class... GridPDEDim1D, class IdxRangeFull, class ExecSpace, class LayoutSpace>
class FFTPoissonSolver<IdxRange<GridPDEDim1D...>, IdxRangeFull, ExecSpace, LayoutSpace>
    : public IPoissonSolver<
              IdxRange<GridPDEDim1D...>,
              IdxRangeFull,
              typename ExecSpace::memory_space,
              LayoutSpace>
{
private:
    using base_type = IPoissonSolver<
            IdxRange<GridPDEDim1D...>,
            IdxRangeFull,
            typename ExecSpace::memory_space,
            LayoutSpace>;

public:
    template <class Dim>
    struct GridFourier : ddc::PeriodicSampling<ddc::Fourier<Dim>>
    {
    };

public:
    using field_type = typename base_type::field_type;
    using const_field_type = typename base_type::const_field_type;

    using vector_field_type = typename base_type::vector_field_type;

    using batch_idx_range_type = typename base_type::batch_idx_range_type;
    using batch_index_type = typename base_type::batch_index_type;

    using laplacian_idx_range_type = typename base_type::laplacian_idx_range_type;

    using layout_space = typename base_type::layout_space;
    using memory_space = typename base_type::memory_space;

    using fourier_idx_range_type
            = IdxRange<GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>;
    using fourier_index_type = typename fourier_idx_range_type::discrete_element_type;

    using fourier_field_mem_type
            = FieldMem<Kokkos::complex<double>, fourier_idx_range_type, memory_space>;
    using fourier_field_type = typename fourier_field_mem_type::span_type;

private:
    static constexpr ddc::FFT_Normalization m_norm = ddc::FFT_Normalization::BACKWARD;

private:
    template <class... FDim>
    KOKKOS_FUNCTION static double get_laplace_operator(Idx<FDim...> index)
    {
        return (((double)ddc::coordinate(ddc::select<FDim>(index))
                 * (double)ddc::coordinate(ddc::select<FDim>(index)))
                + ...);
    }

    template <class Dim>
    void differentiate_and_invert_fourier_values(
            DField<laplacian_idx_range_type, memory_space, LayoutSpace> derivative,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        negative_differentiate_equation<Dim>(fourier_derivative, values);
        // Perform the inverse 1D FFT of the solution to deduce the electric field
        ddc::
                ifft(ExecSpace(),
                     get_field(derivative),
                     get_field(fourier_derivative),
                     ddc::kwArgs_fft {m_norm});
    }

    void get_gradient(
            DField<laplacian_idx_range_type, memory_space, LayoutSpace> gradient,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        using Dim =
                typename ddc::type_seq_element_t<0, ddc::to_type_seq_t<laplacian_idx_range_type>>::
                        continuous_dimension_type;
        using FourierDim = GridFourier<Dim>;
        differentiate_and_invert_fourier_values<FourierDim>(gradient, fourier_derivative, values);
    }

    template <class... Dims>
    void get_gradient(
            VectorField<
                    double,
                    laplacian_idx_range_type,
                    VectorIndexSet<Dims...>,
                    memory_space,
                    layout_space> gradient,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        ((differentiate_and_invert_fourier_values<
                 GridFourier<Dims>>(ddcHelper::get<Dims>(gradient), fourier_derivative, values)),
         ...);
    }

    template <class Grid1D>
    void init_fourier_space(IdxRange<Grid1D> idx_range)
    {
        using GridFFT = GridFourier<typename Grid1D::continuous_dimension_type>;
        ddc::init_discrete_space<GridFFT>(ddc::init_fourier_space<GridFFT>(idx_range));
    }

public:
    template <class Layout>
    void solve_poisson_equation(
            fourier_field_type intermediate_chunk,
            DField<laplacian_idx_range_type, memory_space, Layout> rho) const
    {
        // Compute FFT(rho)
        ddc::fft(ExecSpace(), intermediate_chunk, rho, ddc::kwArgs_fft {m_norm});

        fourier_idx_range_type const k_mesh = get_idx_range(intermediate_chunk);

        // Solve Poisson's equation -\Delta phi = -(\sum_j \partial_j^2) \phi = rho
        //   in Fourier space as -(\sum_j i*k_i * i*k_i) FFT(Phi) = FFT(rho))
        ddc::parallel_for_each(
                ExecSpace(),
                k_mesh,
                KOKKOS_LAMBDA(fourier_index_type const ik) {
                    if (ik != k_mesh.front()) {
                        intermediate_chunk(ik) = intermediate_chunk(ik) / get_laplace_operator(ik);
                    } else {
                        intermediate_chunk(ik) = 0.;
                    }
                });
    }

    template <class Dim>
    void negative_differentiate_equation(fourier_field_type derivative, fourier_field_type values)
            const
    {
        Kokkos::complex<double> imaginary_unit(0.0, 1.0);
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(values),
                KOKKOS_LAMBDA(fourier_index_type const ik) {
                    Idx<Dim> const ikx = ddc::select<Dim>(ik);
                    derivative(ik) = -imaginary_unit * ddc::coordinate(ikx) * values(ik);
                });
    }

public:
    explicit FFTPoissonSolver(laplacian_idx_range_type laplacian_idx_range)
    {
        ((init_fourier_space<GridPDEDim1D>(ddc::select<GridPDEDim1D>(laplacian_idx_range))), ...);
    }

    virtual field_type operator()(field_type phi, field_type rho) const final
    {
        Kokkos::Profiling::pushRegion("FFTPoissonSolver");

        laplacian_idx_range_type idx_range(get_idx_range(phi));
        batch_idx_range_type batch_idx_range(get_idx_range(phi));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);
        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         get_field(intermediate_chunk),
                         ddc::kwArgs_fft {m_norm});
        });

        Kokkos::Profiling::popRegion();
        return phi;
    }

    virtual field_type operator()(field_type phi, vector_field_type E, field_type rho) const final
    {
        Kokkos::Profiling::pushRegion("FFTPoissonSolver");

        laplacian_idx_range_type idx_range(get_idx_range(phi));
        batch_idx_range_type batch_idx_range(get_idx_range(phi));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);
        fourier_field_mem_type fourier_efield_alloc(k_mesh);

        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);
        fourier_field_type fourier_efield = get_field(fourier_efield_alloc);

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);
            get_gradient(E[ib], fourier_efield, intermediate_chunk);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         get_field(intermediate_chunk),
                         ddc::kwArgs_fft {m_norm});
        });
        Kokkos::Profiling::popRegion();
        return phi;
    }
};
```


