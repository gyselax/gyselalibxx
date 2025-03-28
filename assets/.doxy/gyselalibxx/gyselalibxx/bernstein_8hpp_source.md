

# File bernstein.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**bernstein.hpp**](bernstein_8hpp.md)

[Go to the documentation of this file](bernstein_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "cartesian_to_barycentric.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "math_tools.hpp"
#include "view.hpp"

template <class X, class Y, class Corner1Tag, class Corner2Tag, class Corner3Tag, std::size_t D>
class TriangularBernsteinPolynomialBasis
{
public:
    using discrete_dimension_type = TriangularBernsteinPolynomialBasis;

    static constexpr std::size_t rank()
    {
        return 2;
    }

    static constexpr std::size_t degree() noexcept
    {
        return D;
    }

    static constexpr std::size_t nbasis()
    {
        return (D + 1) * (D + 2) / 2;
    }

    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    private:
        CartesianToBarycentric<X, Y, Corner1Tag, Corner2Tag, Corner3Tag> m_coord_changer;

    public:
        using discrete_dimension_type = TriangularBernsteinPolynomialBasis;

        using discrete_element_type = Idx<DDim>;

        using discrete_domain_type = IdxRange<DDim>;

        using discrete_vector_type = IdxStep<DDim>;

        explicit Impl(CartesianToBarycentric<X, Y, Corner1Tag, Corner2Tag, Corner3Tag> const&
                              coord_changer)
            : m_coord_changer(coord_changer)
        {
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_coord_changer(impl.m_coord_changer)
        {
        }

        Impl(Impl const& x) = default;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        Impl& operator=(Impl const& x) = default;

        Impl& operator=(Impl&& x) = default;

        void eval_basis(host_t<DField<IdxRange<DDim>>> values, Coord<X, Y> const& x) const;
    };
};

template <class X, class Y, class Corner1Tag, class Corner2Tag, class Corner3Tag, std::size_t D>
template <class DDim, class MemorySpace>
void TriangularBernsteinPolynomialBasis<X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D>::Impl<
        DDim,
        MemorySpace>::eval_basis(host_t<DField<IdxRange<DDim>>> values, Coord<X, Y> const& x) const
{
    const Coord<Corner1Tag, Corner2Tag, Corner3Tag> bary_coords = m_coord_changer(x);
    const double l1 = ddc::get<Corner1Tag>(bary_coords);
    const double l2 = ddc::get<Corner2Tag>(bary_coords);
    const double l3 = ddc::get<Corner3Tag>(bary_coords);
    assert(values.size() == nbasis());

    Idx<DDim> idx(0);
    for (std::size_t i(0); i < D + 1; ++i) {
        for (std::size_t j(0); j < D + 1 - i; ++j, ++idx) {
            const int k = D - i - j;
            const double multinomial_coefficient
                    = factorial(D) / (factorial(i) * factorial(j) * factorial(k));
            values(idx) = multinomial_coefficient * ipow(l1, i) * ipow(l2, j) * ipow(l3, k);
        }
    }
}
```


