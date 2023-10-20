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
 */
class CollisionsInter : public IRightHandSide
{
private:
    double m_nustar0;
    DFieldSpX m_nustar_profile;

public:
    /**
     * @brief The constructor for the operator.
     *
     * @param[in] mesh The domain on which the operator will act.
     * @param[in] nustar0 The collision coefficient.
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
     * @return The collision coefficient.
     */
    double get_nustar0() const;

private:
    void compute_rhs(DSpanSpXVx rhs, DViewSpXVx allfdistribu) const;
};
