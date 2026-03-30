

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
#include "math_tools.hpp"
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

            int node = check_if_node(x, poly_start);
            if (node == -1) {
                std::array<DataType, D + 1> weights = get_weights(poly_start);
                calculate_values_between_nodes(values, x, poly_start, weights);
            } else {
                for (int j(0); j < D + 1; ++j) {
                    values(j) = static_cast<int>(j == node);
                }
            }
        }

        KOKKOS_INLINE_FUNCTION void eval_basis_and_n_derivs(
                Span2D<DataType> derivs,
                coord_type const& x,
                Idx<knot_grid> poly_start,
                std::size_t n_derivs) const
        {
            KOKKOS_ASSERT(x >= ddc::coordinate(poly_start));
            KOKKOS_ASSERT(x <= ddc::coordinate(poly_start + degree()));
            KOKKOS_ASSERT(n_derivs <= degree());
            KOKKOS_ASSERT(derivs.extent(0) == degree() + 1);
            KOKKOS_ASSERT(derivs.extent(1) == n_derivs + 1);

            constexpr std::size_t n_basis = degree() + 1;
            int node = check_if_node(x, poly_start);
            std::array<DataType, D + 1> weights = get_weights(poly_start);

            // If coordinate not found at a node
            if (node == -1) {
                std::array<DataType, n_basis> vals_ptr;
                Span1D<DataType> values(vals_ptr.data(), n_basis);
                calculate_values_between_nodes(values, x, poly_start, weights);

                std::array<int, n_basis> combinations;

                // For each basis element
                for (std::size_t j(0); j < n_basis; ++j) {
                    derivs(j, 0) = values(j);
                    for (std::size_t n(1); n < n_derivs + 1; ++n) {
                        DataType factor = 0;
                        combinations[0] = -1;
                        for (std::size_t i(1); i < n; ++i) {
                            combinations[i] = 0;
                        }
                        int i = 0;
                        while (i > -1) {
                            combinations[i] += 1;
                            combinations[i] += static_cast<int>(combinations[i] == j);
                            for (std::size_t k(i + 1); k < n; ++k) {
                                combinations[k] = combinations[k - 1] + 1;
                                combinations[k] += static_cast<int>(combinations[k] == j);
                            }
                            i = n - 1;
                            DataType divisor(1);
                            for (std::size_t k(0); k < n; ++k) {
                                divisor *= (x - ddc::coordinate(poly_start + combinations[k]));
                            }
                            factor += 1.0 / divisor;
                            int max_val_along_dim = i + n_basis - n;
                            max_val_along_dim -= static_cast<int>(max_val_along_dim <= j);
                            while (i > -1 and combinations[i] == max_val_along_dim) {
                                i--;
                                max_val_along_dim = i + n_basis - n;
                                max_val_along_dim -= static_cast<int>(max_val_along_dim <= j);
                            }
                        }
                        derivs(j, n) = values(j) * factor * factorial(n);
                    }
                }
            } else {
                DataType xi = ddc::coordinate(poly_start + node);
                for (std::size_t j(0); j < n_basis; ++j) {
                    derivs(j, 0) = static_cast<int>(j == node);
                    DataType xj = ddc::coordinate(poly_start + j);
                    if (j == node)
                        continue;
                    for (std::size_t n(1); n < n_derivs + 1; ++n) {
                        // \partial^nl_j(x_i) = n!/wi [ (-1)^(n+1) wj/(xi-xj)^n
                        //                      - \sum_k=1^{n-1}\sum_{p\ne i} (-1)^(n+1-k)
                        //                                       wp \partial^k l_j(x_i) / k! / (xi-xp)^n ]
                        derivs(j, n) = neg_1_pow(n + 1) * weights[j] / Kokkos::pow(xi - xj, n);
                        for (std::size_t p(0); p < n_basis; ++p) {
                            if (p == node)
                                continue;
                            DataType xp = ddc::coordinate(poly_start + p);
                            for (std::size_t k(1); k < n; ++k) {
                                derivs(j, n) -= neg_1_pow(n + 1 - k) * weights[p] * derivs(j, k)
                                                / factorial(k) / ipow(xi - xp, n - k);
                            }
                        }
                        derivs(j, n) *= factorial(n) / weights[node];
                    }
                }
                for (std::size_t n(1); n < n_derivs + 1; ++n) {
                    derivs(node, n) = 0;
                    for (std::size_t j(0); j < n_basis; ++j) {
                        if (j == node)
                            continue;
                        derivs(node, n) -= derivs(j, n);
                    }
                }
            }
        }

    private:
        KOKKOS_INLINE_FUNCTION static int check_if_node(
                coord_type const& x,
                Idx<knot_grid> poly_start)
        {
            DataType eps = Kokkos::Experimental::epsilon_v<DataType> * 4;
            for (std::size_t i(0); i < D + 1; ++i) {
                if (Kokkos::fabs(x - ddc::coordinate(poly_start + i)) < eps) {
                    return i;
                }
            }
            return -1;
        }

        KOKKOS_INLINE_FUNCTION static std::array<DataType, D + 1> get_weights(
                Idx<knot_grid> poly_start)
        {
            std::array<DataType, D + 1> weights;
            for (std::size_t i(0); i < D + 1; ++i) {
                DataType xi = ddc::coordinate(poly_start + i);
                DataType numerator = 1;
                for (std::size_t j(0); j < D + 1; ++j) {
                    if (i == j)
                        continue;
                    DataType xj = ddc::coordinate(poly_start + j);
                    numerator *= (xi - xj);
                }
                weights[i] = 1.0 / numerator;
            }
            return weights;
        }

        KOKKOS_INLINE_FUNCTION static void calculate_values_between_nodes(
                Span1D<DataType> values,
                coord_type const& x,
                Idx<knot_grid> poly_start,
                std::array<DataType, D + 1> const& weights)
        {
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

template <class DDim>
constexpr bool is_non_uniform_lagrange_basis_v = is_non_uniform_lagrange_basis<DDim>::value;
```


