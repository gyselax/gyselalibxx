// SPDX-License-Identifier: MIT
#pragma once

#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


/**
 * @brief A class which implements a partial derivative operator
 * using a finite differences calculation of order two. A decentered
 * scheme is used at the boundary, whereas centred finite difference
 * are used inside the domain.
 *
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivativeWithBValue
    : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    /// The type of the field to be differentiated.
    using DFieldMemType = DFieldMem<IdxRangeFull>;

    /// The type of a reference to the field to be differentiated.
    using typename base_type::DFieldType;

    /// The type of a constant reference to the field to be differentiated.
    using typename base_type::DConstFieldType;

    /// The index range of the dimension on which the partial derivative is calculated.
    using typename base_type::IdxRangeDeriv;

    /// The index range of all dimensions except DerivativeDimension.
    using typename base_type::IdxRangeBatch;

private:
    /// The field to be differentiated, with the two boundary values
    DFieldMemType m_field;
    double m_b_value_left;
    double m_b_value_right;

public:
    /**
     * @brief Construct an instance of the class CentralFDMPartialDerivative.
     *
     * @param field_ref The field to be differentiated.
     * @param bvalue_left The left boundary value.
     * @param bvalue_right The right boundary value.
     */
    explicit CentralFDMPartialDerivativeWithBValue(
            DConstFieldType const field_ref,
            double bvalue_left,
            double bvalue_right)
        : m_field(get_idx_range(field_ref))
        , m_b_value_left(bvalue_left)
        , m_b_value_right(bvalue_right)
    {
        ddc::parallel_deepcopy(get_field(m_field), field_ref);
    }

    /**
     * @brief Compute the partial derivative of a field in a given direction
     * using a finite difference scheme. For more information about the coefficients,
     * see `./README.md`
     *
     * @param[out] differentiated_field On output, contains values of the differentiated field.
     */
    void operator()(DFieldType differentiated_field) const final
    {
        using IdxFull = typename IdxRangeFull::discrete_element_type;
        using IdxDeriv = typename IdxRangeDeriv::discrete_element_type;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;
        using IdxStepDeriv = typename IdxRangeDeriv::discrete_vector_type;

        IdxRangeFull idxrange_full = get_idx_range(m_field);
        IdxRangeDeriv idxrange_deriv(idxrange_full);
        IdxRangeBatch idxrange_batch(idxrange_full);

        DConstFieldType const field_proxy = get_const_field(m_field);

        IdxStepDeriv const step(1);

        // front batched derivative
        // we compute second order left and right decentred FDM (using boundary values)
        // and we take the average
        {
            IdxDeriv const ix(idxrange_deriv.front());
            double const h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
            double const h2 = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
            double const c3 = -h1 / (h2 * (h1 + h2));
            double const c2 = 1. / h1 + 1. / h2;
            double bvalue_left(m_b_value_left);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        double value_left = c2 * field_proxy(ib, ix + step)
                                            + c3 * field_proxy(ib, ix + 2 * step);
                        double value_right = -c2 * bvalue_left - c3 * bvalue_left;
                        differentiated_field(ibx) = (value_right + value_left) / 2;
                    });
        }

        // back batched derivative
        {
            IdxDeriv const ix = idxrange_deriv.back();
            double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
            double const h2 = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
            double const c3 = h1 / (h2 * (h1 + h2));
            double const c2 = -(h1 + h2) / (h1 * h2);
            double bvalue_right(m_b_value_right);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        double value_left = c2 * field_proxy(ib, ix - step)
                                            + c3 * field_proxy(ib, ix - 2 * step);
                        double value_right = -c2 * bvalue_right - c3 * bvalue_right;
                        differentiated_field(ibx) = (value_left + value_right) / 2;
                    });
        }

        // central domain batched derivative
        {
            IdxRangeDeriv idxrange_deriv_central = idxrange_deriv.remove(step, step);
            IdxRangeFull idxrange_central(idxrange_deriv_central, idxrange_batch);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_central,
                    KOKKOS_LAMBDA(IdxFull ibx) {
                        IdxBatch ib(ibx);
                        IdxDeriv ix(ibx);
                        double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                        double const h2 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                        double const c3 = h1 / (h2 * (h1 + h2));
                        double const c2 = 1. / h1 - 1. / h2;
                        double const c1 = -c3 - c2;
                        differentiated_field(ibx) = c1 * field_proxy(ib, ix - step)
                                                    + c2 * field_proxy(ibx)
                                                    + c3 * field_proxy(ib, ix + step);
                    });
        }
    }
};

/**
 * @brief A class which stores information necessary to create a pointer to 
 * an instance of the CentralFDMPartialDerivativeWithBValue class.
 *
 * This class allows an instance of the CentralFDMPartialDerivativeWithBValue class to be instantiated where necessary.
 * Typically, the CentralFDMPartialDerivativeCreator is instantiated in the initialisation of the simulation, 
 * and the corresponding CentralFDMPartialDerivativeWithBValue object is instantiated where computing partial derivatives
 * is required. 
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivativeWithBValueCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    /// The type of a constant reference to the field to be differentiated.
    using DConstFieldType = DConstField<IdxRangeFull>;

private:
    // The boundary values for the derivative
    double m_b_value_left;
    double m_b_value_right;

public:
    /**
     * @brief Construct an instance of the CentralFDMPartialDerivativeWithBValueCreator class.
     *
     * @param bvalue_left The left boundary value.
     * @param bvalue_right The right boundary value.
     */
    CentralFDMPartialDerivativeWithBValueCreator(double bvalue_left, double bvalue_right)
        : m_b_value_left(bvalue_left)
        , m_b_value_right(bvalue_right)
    {
    }

    /**
     * Create a pointer to an instance of the abstract class IPartialDerivative.
     * The type of the returned object will be determined when the pointer is 
     * dereferenced.
     *
     * @param[in] field_ref A field to be differentiated.
     *
     * @return A pointer to an instance of the IPartialDerivative class.
     */
    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstFieldType field_ref) const final
    {
        return std::make_unique<CentralFDMPartialDerivativeWithBValue<
                IdxRangeFull,
                DerivativeDimension>>(field_ref, m_b_value_left, m_b_value_right);
    }
};
