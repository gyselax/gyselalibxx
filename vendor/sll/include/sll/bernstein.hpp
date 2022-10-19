#pragma once
#include <ddc/ddc.hpp>

#include <sll/mapping/barycentric_coordinates.hpp>
#include <sll/math_tools.hpp>
#include <sll/view.hpp>

template <
        class Tag1,
        class Tag2,
        class Corner1Tag,
        class Corner2Tag,
        class Corner3Tag,
        std::size_t D>
class BernsteinPolynomialBasis
{
public:
    using discrete_element_type = DiscreteElement<BernsteinPolynomialBasis>;

    using discrete_domain_type = DiscreteDomain<BernsteinPolynomialBasis>;

    using discrete_vector_type = DiscreteVector<BernsteinPolynomialBasis>;

public:
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

    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        CartesianToBarycentricCoordinates<Tag1, Tag2, Corner1Tag, Corner2Tag, Corner3Tag>
                m_coord_changer;

    public:
        using discrete_dimension_type = BernsteinPolynomialBasis;

        Impl(CartesianToBarycentricCoordinates<
                Tag1,
                Tag2,
                Corner1Tag,
                Corner2Tag,
                Corner3Tag> const& coord_changer)
            : m_coord_changer(coord_changer)
        {
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<OriginMemorySpace> const& impl) : m_coord_changer(impl.m_coord_changer)
        {
        }

        Impl(Impl const& x) = default;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        Impl& operator=(Impl const& x) = default;

        Impl& operator=(Impl&& x) = default;

        void eval_basis(
                ChunkSpan<double, DiscreteDomain<BernsteinPolynomialBasis>> values,
                Coordinate<Tag1, Tag2> const& x) const;
    };
};

template <
        class Tag1,
        class Tag2,
        class Corner1Tag,
        class Corner2Tag,
        class Corner3Tag,
        std::size_t D>
template <class MemorySpace>
void BernsteinPolynomialBasis<Tag1, Tag2, Corner1Tag, Corner2Tag, Corner3Tag, D>::
        Impl<MemorySpace>::eval_basis(
                ChunkSpan<double, DiscreteDomain<BernsteinPolynomialBasis>> values,
                Coordinate<Tag1, Tag2> const& x) const
{
    const Coordinate<Corner1Tag, Corner2Tag, Corner3Tag> bary_coords = m_coord_changer(x);
    const double l1 = get<Corner1Tag>(bary_coords);
    const double l2 = get<Corner2Tag>(bary_coords);
    const double l3 = get<Corner3Tag>(bary_coords);
    assert(values.size() == nbasis());

    DiscreteElement<BernsteinPolynomialBasis> idx(0);
    for (std::size_t i(0); i < D + 1; ++i) {
        for (std::size_t j(0); j < D + 1 - i; ++j, ++idx) {
            const int k = D - i - j;
            const double multinomial_coefficient
                    = factorial(D) / (factorial(i) * factorial(j) * factorial(k));
            values(idx) = multinomial_coefficient * ipow(l1, i) * ipow(l2, j) * ipow(l3, k);
        }
    }
}
