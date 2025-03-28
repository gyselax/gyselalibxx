

# File inv\_jacobian\_o\_point.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**inv\_jacobian\_o\_point.hpp**](inv__jacobian__o__point_8hpp.md)

[Go to the documentation of this file](inv__jacobian__o__point_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_to_cartesian.hpp"
#include "indexed_tensor.hpp"
#include "mapping_tools.hpp"
#include "view.hpp"

template <class Mapping, class CoordRTheta>
class InvJacobianOPoint;

// Pre-declaration of CombinedMapping.
template <class Mapping1, class Mapping2>
class CombinedMapping;

template <class X, class Y, class R, class Theta, class Xpc, class Ypc>
class InvJacobianOPoint<
        CombinedMapping<
                CircularToCartesian<R, Theta, X, Y>,
                CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    using CoordRTheta = Coord<R, Theta>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;

private:
    using Mapping = CombinedMapping<
            CircularToCartesian<R, Theta, X, Y>,
            CartesianToCircular<Xpc, Ypc, R, Theta>>;

private:
    Mapping m_mapping;

public:
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>>
    operator()() const
    {
        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J;
        ddcHelper::get<Xpc, X_cov>(J) = 1.;
        ddcHelper::get<Xpc, Y_cov>(J) = 0.;
        ddcHelper::get<Ypc, X_cov>(J) = 0.;
        ddcHelper::get<Ypc, Y_cov>(J) = 1.;
        return J;
    }
};

template <class X, class Y, class R, class Theta, class Xpc, class Ypc>
class InvJacobianOPoint<
        CombinedMapping<CzarnyToCartesian<R, Theta, X, Y>, CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    using CoordRTheta = Coord<R, Theta>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;

private:
    using CzarnyToCart = CzarnyToCartesian<R, Theta, X, Y>;
    using Mapping = CombinedMapping<CzarnyToCart, CartesianToCircular<Xpc, Ypc, R, Theta>>;

private:
    Mapping m_mapping;

public:
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>>
    operator()() const
    {
        const double epsilon = m_mapping.template get<CzarnyToCart>().epsilon();
        const double e = m_mapping.template get<CzarnyToCart>().e();
        const double xi = Kokkos::sqrt(1. / (1. - epsilon * epsilon * 0.25));
        const double sqrt_eps_2 = Kokkos::sqrt(1. + epsilon * epsilon);

        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J;
        ddcHelper::get<Xpc, X_cov>(J) = -sqrt_eps_2;
        ddcHelper::get<Xpc, Y_cov>(J) = 0.;
        ddcHelper::get<Ypc, X_cov>(J) = 0.;
        ddcHelper::get<Ypc, Y_cov>(J) = (2 - sqrt_eps_2) / e / xi;
        return J;
    }
};

template <
        class X,
        class Y,
        class SplineEvaluator,
        class R,
        class Theta,
        class MemorySpace,
        class Xpc,
        class Ypc>
class InvJacobianOPoint<
        CombinedMapping<
                DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>,
                CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    using CoordRTheta = Coord<R, Theta>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using Xpc_cov = typename Xpc::Dual;
    using Ypc_cov = typename Ypc::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    using Mapping = CombinedMapping<
            DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>,
            CartesianToCircular<Xpc, Ypc, R, Theta>>;

    using IdxRangeR = typename SplineEvaluator::evaluation_domain_type1;
    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxRangeRTheta = typename SplineEvaluator::evaluation_domain_type;
    using IdxR = typename IdxRangeR::discrete_element_type;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    Mapping m_mapping;

public:
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    KOKKOS_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> operator()()
            const
    {
        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace> const& discrete_mapping
                = m_mapping.template get<
                        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>();
        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J(0.0);
        IdxRangeRTheta idx_range_singular_point = discrete_mapping.idx_range_singular_point();
        // Average the values at (r = 0, theta):
        IdxR ir(idx_range_singular_point.front());
        for (IdxTheta itheta : IdxRangeTheta(idx_range_singular_point)) {
            Coord<R, Theta> coord(ddc::coordinate(ir), ddc::coordinate(itheta));
            DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> J_first_order
                    = discrete_mapping.first_order_jacobian_matrix_r_rtheta(coord);

            double th = ddc::get<Theta>(coord);
            DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<Xpc_cov, Ypc_cov>> J_circ_r_rtheta;
            ddcHelper::get<R, Xpc_cov>(J_circ_r_rtheta) = Kokkos::cos(th);
            ddcHelper::get<R, Ypc_cov>(J_circ_r_rtheta) = Kokkos::sin(th);
            ddcHelper::get<Theta, Xpc_cov>(J_circ_r_rtheta) = -Kokkos::sin(th);
            ddcHelper::get<Theta, Ypc_cov>(J_circ_r_rtheta) = Kokkos::cos(th);

            Tensor J_theta = inverse(
                    tensor_mul(index<'i', 'j'>(J_first_order), index<'j', 'k'>(J_circ_r_rtheta)));

            J += J_theta;
        }

        int const theta_size = idx_range_singular_point.size();
        J /= theta_size;
        return J;
    }
};
```


