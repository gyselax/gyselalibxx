

# File lagrange\_basis\_non\_uniform.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**lagrange\_basis\_non\_uniform.hpp**](lagrange__basis__non__uniform_8hpp.md)

[Go to the documentation of this file](lagrange__basis__non__uniform_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "ddc_aliases.hpp"
#include "view.hpp"

namespace detail {

class NonUniformLagrangeBasisBase
{
};

} // namespace detail

template <class T>
struct NonUniformLagrangeKnots : NonUniformGridBase<typename T::continuous_dimension_type>
{
};

template <class Dim, std::size_t D, class DataType = double>
class NonUniformLagrangeBasis : detail::NonUniformLagrangeBasisBase
{
    static_assert(D > 0, "Parameter `D` must be positive");
    static_assert(std::is_floating_point_v<DataType>);

public:
    using continuous_dimension_type = Dim;

    using coord_type = Coord<Dim>;

    using discrete_dimension_type = NonUniformLagrangeBasis;

    static constexpr std::size_t degree() noexcept
    {
        return D;
    }

    static constexpr bool is_periodic() noexcept
    {
        return continuous_dimension_type::PERIODIC;
    }

    static constexpr bool is_uniform() noexcept
    {
        return false;
    }


    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    public:
        using knot_grid = NonUniformLagrangeKnots<DDim>;

    private:
        IdxRange<knot_grid> m_knot_domain;
        IdxRange<knot_grid> m_break_point_domain;

        Idx<DDim> m_reference;

    public:
        Impl() = default;

        template <class Grid1D>
        explicit Impl(IdxRange<Grid1D> break_point_domain)
            : m_reference(ddc::create_reference_discrete_element<DDim>())
        {
            static_assert(std::is_same_v<typename Grid1D::continuous_dimension_type, Dim>);
            static_assert(
                    ddc::is_uniform_point_sampling<Grid1D>()
                    || ddc::is_non_uniform_point_sampling<Grid1D>());
            assert(ddc::is_discrete_space_initialized<Grid1D>());
            assert(break_point_domain.size() >= D + 1);

            // Initialise knot grid
            if constexpr (is_periodic()) {
                std::size_t npoints = ddc::discrete_space<Grid1D>().size();
                std::vector<coord_type> points(npoints + D);
                Idx<Grid1D> idx_front = ddc::discrete_space<Grid1D>().front();
                Idx<Grid1D> idx_back = idx_front + npoints - 1;
                for (std::size_t i(0); i < npoints; ++i) {
                    points[i] = ddc::coordinate(idx_front + i);
                }
                for (std::size_t i(1); i < D + 1; ++i) {
                    points[npoints + i - 1]
                            = ddc::coordinate(idx_back)
                              + (ddc::coordinate(idx_front + i) - ddc::coordinate(idx_front));
                }
                ddc::init_discrete_space<knot_grid>(points);
                m_knot_domain = IdxRange<knot_grid>(
                        ddc::discrete_space<knot_grid>().front(),
                        IdxStep<knot_grid>(break_point_domain.size() + D));
                m_break_point_domain = m_knot_domain.remove_last(IdxStep<knot_grid>(D));
            } else {
                std::size_t npoints = ddc::discrete_space<Grid1D>().size();
                Idx<Grid1D> idx_front = ddc::discrete_space<Grid1D>().front();
                std::vector<coord_type> points(npoints);
                for (std::size_t i(0); i < npoints; ++i) {
                    points[i] = ddc::coordinate(idx_front + i);
                }
                ddc::init_discrete_space<knot_grid>(points);
                m_knot_domain = IdxRange<knot_grid>(
                        ddc::discrete_space<knot_grid>().front(),
                        IdxStep<knot_grid>(break_point_domain.size()));
                m_break_point_domain = m_knot_domain;
            }
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_knot_domain(impl.m_knot_domain)
            , m_break_point_domain(impl.m_break_point_domain)
            , m_reference(impl.m_reference)
        {
        }

        KOKKOS_INLINE_FUNCTION coord_type rmin() const noexcept
        {
            return ddc::coordinate(m_break_point_domain.front());
        }

        KOKKOS_INLINE_FUNCTION coord_type rmax() const noexcept
        {
            return ddc::coordinate(m_break_point_domain.back());
        }

        KOKKOS_INLINE_FUNCTION DataType length() const noexcept
        {
            return static_cast<DataType>(rmax()) - static_cast<DataType>(rmin());
        }

        KOKKOS_INLINE_FUNCTION IdxRange<knot_grid> break_point_domain() const
        {
            return m_break_point_domain;
        }

        KOKKOS_INLINE_FUNCTION IdxRange<knot_grid> full_domain() const
        {
            return m_knot_domain;
        }

        KOKKOS_INLINE_FUNCTION void eval_basis(
                Span1D<DataType> values,
                coord_type const& x,
                Idx<knot_grid> poly_start) const
        {
            KOKKOS_ASSERT(values.size() == degree() + 1);
            KOKKOS_ASSERT(x >= ddc::coordinate(poly_start));
            KOKKOS_ASSERT(x <= ddc::coordinate(poly_start + degree()));

            DataType eps = Kokkos::Experimental::epsilon_v<DataType> * 4;
            for (std::size_t i(0); i < D + 1; ++i) {
                if (Kokkos::fabs(x - ddc::coordinate(poly_start + i)) < eps) {
                    for (int j(0); j < D + 1; ++j) {
                        values(j) = static_cast<int>(j == i);
                    }
                    return;
                }
            }
            std::array<DataType, D + 1> weights;
            for (std::size_t i(0); i < D + 1; ++i) {
                DataType numerator = 1;
                for (std::size_t j(0); j < D + 1; ++j) {
                    if (i == j)
                        continue;
                    numerator
                            *= (ddc::coordinate(poly_start + i) - ddc::coordinate(poly_start + j));
                }
                weights[i] = 1.0 / numerator;
            }
            for (int i(0); i < D + 1; ++i) {
                DataType numerator = weights[i] / (x - ddc::coordinate(poly_start + i));
                DataType denominator(0.0);
                for (int j(0); j < D + 1; ++j) {
                    denominator += weights[j] / (x - ddc::coordinate(poly_start + j));
                }
                values(i) = numerator / denominator;
            }
        }
    };
};

template <class DDim>
struct is_non_uniform_lagrange_basis
    : public std::is_base_of<detail::NonUniformLagrangeBasisBase, DDim>::type
{
};
```


