// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "cartesian_to_barycentric.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "math_tools.hpp"
#include "view.hpp"

/**
 * @brief A class which evaluates the triangular Bernstein polynomials.
 *
 * Triangular Bernstein polynomials of degree @f$ D @f$ are defined as:
 * @f$ T_{i_1,i_2,i_3}(x,y) = \binom{D}{i_1\,i_2\,i_3} \lambda_1^{i_1} \lambda_2^{i_2} \lambda_3^{i_3}, \quad i_1 +i_2+i_3 =D @f$
 *
 * Where @f$ \lambda_1 @f$, @f$ \lambda_2 @f$, and @f$ \lambda_3 @f$ are the barycentric coordinates of @f$ (x,y) @f$
 *
 * c.f. Toshniwal et al. 2017
 * Multi-degree smooth polar splines: A framework for geometric modeling and isogeometric analysis
 * https://doi.org/10.1016/j.cma.2016.11.009
 *
 * @tparam X The first dimension of the Cartesian coordinate system.
 * @tparam Y The second dimension of the Cartesian coordinate system.
 * @tparam Corner1Tag A class to identify the first corner.
 * @tparam Corner2Tag A class to identify the second corner.
 * @tparam Corner3Tag A class to identify the third corner.
 * @tparam D The degree of the polynomial.
 */
template <class X, class Y, class Corner1Tag, class Corner2Tag, class Corner3Tag, std::size_t D>
class TriangularBernsteinPolynomialBasis
{
public:
    /// A tag for DDC to recognise the discrete dimension type.
    using discrete_dimension_type = TriangularBernsteinPolynomialBasis;

    /**
     * @brief The rank of the system of equations. This is equal to the number of dimensions
     * in the coordinates on which the polynomials are evaluated.
     *
     * @return The rank of the system.
     */
    static constexpr std::size_t rank()
    {
        return 2;
    }

    /**
     * @brief The degree of the triangular Bernstein polynomials.
     * @return The degree of the triangular Bernstein polynomials.
     */
    static constexpr std::size_t degree() noexcept
    {
        return D;
    }

    /**
     * @brief The number of basis elements.
     * @return The number of basis elements.
     */
    static constexpr std::size_t nbasis()
    {
        return (D + 1) * (D + 2) / 2;
    }

    /**
     * The Impl class holds the implementation of the TriangularBernsteinPolynomialBasis.
     *
     * @tparam MemorySpace Indicates where the object is saved. This is either on the host or the device.
     */
    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    private:
        CartesianToBarycentric<X, Y, Corner1Tag, Corner2Tag, Corner3Tag> m_coord_changer;

    public:
        /// The tag which identifies the basis.
        using discrete_dimension_type = TriangularBernsteinPolynomialBasis;

        /// The type of an index of an element of the basis.
        using discrete_element_type = Idx<DDim>;

        /// The type of the index range of the basis.
        using discrete_domain_type = IdxRange<DDim>;

        /// The type of an index step from one element of the basis to another.
        using discrete_vector_type = IdxStep<DDim>;

        /**
         * @brief Construct the basis from the barycentric coordinate mapping.
         * @param[in] coord_changer The class which converts Cartesian coordinates to
         * barycentric coordinates.
         */
        explicit Impl(CartesianToBarycentric<X, Y, Corner1Tag, Corner2Tag, Corner3Tag> const&
                              coord_changer)
            : m_coord_changer(coord_changer)
        {
        }

        /**
         * @brief Construct the basis by copy. This constructor is used to create
         * the class on a different memory space.
         * @param[in] impl The implementation of the origin memory space.
         */
        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_coord_changer(impl.m_coord_changer)
        {
        }

        /**
         * @brief Construct the basis by copy.
         * @param[in] x The basis to be copied.
         */
        Impl(Impl const& x) = default;

        /**
         * @brief Construct the basis from an r-value.
         * @param[in] x The temporary basis to be copied.
         */
        Impl(Impl&& x) = default;

        ~Impl() = default;

        /**
         * @brief Copy-assign the class.
         *
         * @param x An reference to another Impl.
         * @return A reference to this object.
         */
        Impl& operator=(Impl const& x) = default;

        /**
         * @brief Move-assign the class.
         *
         * @param x An rvalue to another Impl.
         * @return A reference to this object.
         */
        Impl& operator=(Impl&& x) = default;

        /**
         * @brief Evaluate the basis at the given coordinate.
         * @param[out] values The values of the basis functions at the coordinate.
         * @param[in] x The coordinate where the polynomials should be evaluated.
         */
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
