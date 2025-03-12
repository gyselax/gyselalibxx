// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>

#include "ddc/discrete_domain.hpp"
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
class CentralFDMPartialDerivative : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    /// The type of a reference to the field to be differentiated.
    using typename base_type::DFieldType;

    /// The type of a constant reference to the field to be differentiated.
    using typename base_type::DConstFieldType;

    /// The index range of the dimension on which the partial derivative is calculated.
    using typename base_type::IdxRangeDeriv;

    /// The index range of all dimensions except DerivativeDimension.
    using typename base_type::IdxRangeBatch;

    /// A constant reference to the field to be differentiated
    DConstFieldType const m_field;

public:
    /**
     * @brief Construct an instance of the class CentralFDMPartialDerivative.
     *
     * @param field The field to be differentiated.
     */
    explicit CentralFDMPartialDerivative(DConstFieldType const field) : m_field(field) {}

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

        DConstFieldType const field_proxy = m_field;

        // front batched derivative
        IdxDeriv const ix_front(idxrange_deriv.front());
        IdxStepDeriv const step(1);
        double const h1_front = ddc::coordinate(ix_front + step) - ddc::coordinate(ix_front);
        double const h2_front
                = ddc::coordinate(ix_front + 2 * step) - ddc::coordinate(ix_front + step);
        double const c3_front = -h1_front / (h2_front * (h1_front + h2_front));
        double const c2_front = 1. / h1_front + 1. / h2_front;
        double const c1_front = -c3_front - c2_front;
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_batch,
                KOKKOS_LAMBDA(IdxBatch ib) {
                    IdxFull ibx(ib, ix_front);
                    differentiated_field(ibx) = c1_front * field_proxy(ibx)
                                                + c2_front * field_proxy(ib, ix_front + step)
                                                + c3_front * field_proxy(ib, ix_front + 2 * step);
                });

        // back batched derivative
        IdxDeriv const ix_back = idxrange_deriv.back();
        double const h1_back = ddc::coordinate(ix_back) - ddc::coordinate(ix_back - step);
        double const h2_back
                = ddc::coordinate(ix_back - step) - ddc::coordinate(ix_back - 2 * step);
        double const c3_back = h1_back / (h2_back * (h1_back + h2_back));
        double const c2_back = -(h1_back + h2_back) / (h1_back * h2_back);
        double const c1_back = -c3_back - c2_back;
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_batch,
                KOKKOS_LAMBDA(IdxBatch ib) {
                    IdxFull ibx(ib, ix_back);
                    differentiated_field(ibx) = c1_back * field_proxy(ibx)
                                                + c2_back * field_proxy(ib, ix_back - step)
                                                + c3_back * field_proxy(ib, ix_back - 2 * step);
                });

        // central domain batched derivative
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
};

/**
 * @brief A class which stores information necessary to create a pointer to 
 * an instance of the CentralFDMPartialDerivative class.
 *
 * This class allows an instance of the CentralFDMPartialDerivative class to be instantiated where necessary.
 * Typically, the CentralFDMPartialDerivativeCreator is instantiated in the initialisation of the simulation, 
 * and the corresponding CentralFDMPartialDerivative object is instantiated where computing partial derivatives
 * is required. 
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    /// The type of a constant reference to the field to be differentiated.
    using DConstFieldType = DConstField<IdxRangeFull>;

public:
    /**
     * Create a pointer to an instance of the abstract class IPartialDerivative.
     * The type of the returned object will be determined when the pointer is 
     * dereferenced.
     *
     * @param[in] field A field to be differentiated.
     *
     * @return A pointer to an instance of the IPartialDerivative class.
     */
    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstFieldType field) const final
    {
        return std::make_unique<CentralFDMPartialDerivative<IdxRangeFull, DerivativeDimension>>(
                field);
    }
};
