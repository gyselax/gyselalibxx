// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
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

    using IdxFull = typename IdxRangeFull::discrete_element_type;
    using IdxDeriv = typename IdxRangeDeriv::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxStepDeriv = typename IdxRangeDeriv::discrete_vector_type;

public:
    static constexpr int n_local_points()
    {
        return 3;
    }

private:
    /// The field to be differentiated
    DFieldMemType m_field;

public:
    /**
     * @brief Construct an instance of the class CentralFDMPartialDerivative.
     *
     * @param field_ref The field to be differentiated.
     */
    explicit CentralFDMPartialDerivative(DConstFieldType const field_ref)
        : m_field(get_idx_range(field_ref))
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
        IdxRangeFull idxrange_full = get_idx_range(m_field);
        IdxRangeDeriv idxrange_deriv(idxrange_full);
        IdxRangeBatch idxrange_batch(idxrange_full);

        DConstFieldType const field_proxy = get_const_field(m_field);

        IdxStepDeriv const step(1);

        if constexpr (!DerivativeDimension::PERIODIC) {
            // front batched derivative
            {
                IdxDeriv const ix(idxrange_deriv.front());
                double const h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                double const h2 = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
                double const c3 = -h1 / (h2 * (h1 + h2));
                double const c2 = 1. / h1 + 1. / h2;
                double const c1 = -c3 - c2;
                const std::source_location location = std::source_location::current();
                ddc::parallel_for_each(
                        location.function_name(),
                        Kokkos::DefaultExecutionSpace(),
                        idxrange_batch,
                        KOKKOS_LAMBDA(IdxBatch ib) {
                            IdxFull ibx(ib, ix);
                            differentiated_field(ibx) = c1 * field_proxy(ibx)
                                                        + c2 * field_proxy(ib, ix + step)
                                                        + c3 * field_proxy(ib, ix + 2 * step);
                        });
            }

            // back batched derivative
            {
                IdxDeriv const ix = idxrange_deriv.back();
                double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                double const h2 = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
                double const c3 = h1 / (h2 * (h1 + h2));
                double const c2 = -(h1 + h2) / (h1 * h2);
                double const c1 = -c3 - c2;
                const std::source_location location = std::source_location::current();
                ddc::parallel_for_each(
                        location.function_name(),
                        Kokkos::DefaultExecutionSpace(),
                        idxrange_batch,
                        KOKKOS_LAMBDA(IdxBatch ib) {
                            IdxFull ibx(ib, ix);
                            differentiated_field(ibx) = c1 * field_proxy(ibx)
                                                        + c2 * field_proxy(ib, ix - step)
                                                        + c3 * field_proxy(ib, ix - 2 * step);
                        });
            }
        } else {
            IdxDeriv const ix_front(idxrange_deriv.front());
            IdxDeriv const ix_back(idxrange_deriv.back());
            double domain_length = ddcHelper::total_interval_length(idxrange_deriv);
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        differentiated_field(ib, ix_back) = fdm_centred(
                                field_proxy(ib, ix_back - step),
                                field_proxy(ib, ix_back),
                                field_proxy(ib, ix_front),
                                ddc::coordinate(ix_back - step),
                                ddc::coordinate(ix_back),
                                ddc::coordinate(ix_front) + domain_length);
                        differentiated_field(ib, ix_front) = fdm_centred(
                                field_proxy(ib, ix_back),
                                field_proxy(ib, ix_front),
                                field_proxy(ib, ix_front + step),
                                ddc::coordinate(ix_back) - domain_length,
                                ddc::coordinate(ix_front),
                                ddc::coordinate(ix_front + step));
                    });
        }

        // central domain batched derivative
        {
            IdxRangeDeriv idxrange_deriv_central = idxrange_deriv.remove(step, step);
            IdxRangeFull idxrange_central(idxrange_deriv_central, idxrange_batch);
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_central,
                    KOKKOS_LAMBDA(IdxFull ibx) {
                        IdxBatch ib(ibx);
                        IdxDeriv ix(ibx);
                        differentiated_field(ibx) = fdm_centred(
                                field_proxy(ib, ix - step),
                                field_proxy(ib, ix),
                                field_proxy(ib, ix + step),
                                ddc::coordinate(ix - step),
                                ddc::coordinate(ix),
                                ddc::coordinate(ix + step));
                    });
        }
    }

    /**
     * @brief Compute the partial derivative locally from the values at 3 known positions.
     *
     * @param[in] field_positions The positions where the field is provied, in ascending order.
     * @param[in] field_values The values of the field at the provided positions.
     * @return The derivative at the central point.
     */
    static KOKKOS_INLINE_FUNCTION double get_derivative(
            std::array<Coord<DerivativeDimension>, 3> const& field_positions,
            std::array<double, 3> const& field_values)
    {
        return fdm_centred(
                field_values[0],
                field_values[1],
                field_values[2],
                field_positions[0],
                field_positions[1],
                field_positions[2]);
    }

private:
    static KOKKOS_INLINE_FUNCTION double fdm_centred(
            double const field1,
            double const field2,
            double const field3,
            double const x1,
            double const x2,
            double const x3)
    {
        double const h1 = x2 - x1;
        double const h2 = x3 - x2;
        double const c3 = h1 / (h2 * (h1 + h2));
        double const c2 = 1. / h1 - 1. / h2;
        double const c1 = -c3 - c2;
        return c1 * field1 + c2 * field2 + c3 * field3;
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
    /// The type of the partial derivative calculator that will be produced.
    using partial_derivative_type = CentralFDMPartialDerivative<IdxRangeFull, DerivativeDimension>;

public:
    /**
     * @brief Construct an instance of the CentralFDMPartialDerivativeCreator class.
     */
    CentralFDMPartialDerivativeCreator() = default;

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
        return std::make_unique<CentralFDMPartialDerivative<IdxRangeFull, DerivativeDimension>>(
                field_ref);
    }
};
