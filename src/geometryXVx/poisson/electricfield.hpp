// SPDX-License-Identifier: MIT

#include <geometry.hpp>
#include <species_info.hpp>


//===========================================================================
// Compute efield = -dPhi/dx where Phi is the electrostatic potential
//  input : Phi values
//===========================================================================
/**
 * @brief An operator which computes the electric field using splines derivation.
 *
 * An operator which computes the electric field:
 * @f$ E = - \frac{d\phi}{dx} @f$ where Phi is the electrostatic potential.
 * This operator uses spline derivation.
 */
class ElectricField
{
    SplineXBuilder_1d const& m_spline_x_builder;

    SplineXEvaluator_1d m_spline_x_evaluator;

public:
    /**
     * Construct the ElectricField operator.
     *
     * @param spline_x_builder A spline builder which calculates the coefficients of a spline representation.
     * @param spline_x_evaluator A spline evaluator which provides the value of a spline representation from its coefficients.
     */
    ElectricField(
            SplineXBuilder_1d const& spline_x_builder,
            SplineXEvaluator_1d const& spline_x_evaluator);

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electric_field The electric_field, the result of the operation.
     * @param[out] electrostatic_potential The electrostatic potential, the input of the operator.
     */
    void operator()(DSpanX electric_field, DViewX electrostatic_potential) const;

    /**
     * The operator which solves the equation using the method described by the class.
     *
     * @param[out] electric_field The electric_field, the result of the operation.
     * @param[out] electrostatic_potential The electrostatic potential, the input of the operator.
     */
    void operator()(DSpanX electric_field, DBSViewX electrostatic_potential) const;
};
