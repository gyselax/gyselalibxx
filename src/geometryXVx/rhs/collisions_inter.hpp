#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

/**
 * @brief Class describing the inter-species collision operator
 * 
 * The inter-species collision operator accounts for momentum and 
 * energy transfer between the maxwellian parts of the distribution
 * function of different species. It is solved using a explicit time 
 * integrator (RK2 for instance).
 * 
 * The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/src/geometryXVx/rhs/doc/collisions_intra_inter.pdf). 
 */
class CollisionsInter : public IRightHandSide
{
private:
    double m_nustar0;
    DFieldSpX m_nustar_profile_alloc;
    DSpanSpX m_nustar_profile;

public:
    /**
     * @brief The constructor for the operator.
     *
     * @param[in] mesh The domain on which the operator will act.
     * @param[in] nustar0 The normalized collisionality.
     */
    CollisionsInter(IDomainSpXVx const& mesh, double nustar0);

    ~CollisionsInter() = default;

    /**
     * @brief Update the distribution function for inter-species collision.
     *
     * Update the distribution function for both electrons and ions to show how
     * it is modified following collisions between the various species.
     * This operator only handles collisions between particles of different
     * species.
     *
     * @param[inout] allfdistribu The distribution function.
     * @param[in] dt The time step over which the collisions occur.
     *
     * @return A span referencing the distribution function passed as argument.
     */
    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

    /**
     * @brief Get the collision coefficient.
     *
     * @return The collisionality.
     */
    double get_nustar0() const;

    /**
     * @brief Computes the expression of the time derivative of the distribution function.
     *  
     * The expression is df = C, where C is the inter species collision operator. 
     * This function is made for the time integrator that is used to solve the 
     * collision operator (RK2 for instance).
     *
     * @param[inout] df The time derivative.
     * @param[in] allfdistribu The distribution function.
     */
    void get_derivative(DSpanSpXVx df, DViewSpXVx allfdistribu) const;
};
