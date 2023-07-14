#pragma once

#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include "geometry.hpp"
#include "i_interpolator_2d_rp.hpp"


/**
 * @brief A class for interpolating a function using splines in polar coordinates.
 *
 * @tag Spline_interpolator_polar
 */
class SplineInterpolatorRP : public IInterpolatorRP
{
private:
    SplineRPBuilder const& m_builder;

    SplineRPEvaluator const& m_evaluator;

    mutable ddc::Chunk<double, BSDomainRP> m_coefs;

public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolatorRP(SplineRPBuilder const& builder, SplineRPEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.spline_domain())
    {
    }

    ~SplineInterpolatorRP() override = default;

    /**
     * @brief Approximate the value of a function at a set of polar coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data
     * 					On input: an array containing the value of the function at the interpolation points.
     * 					On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates
     * 					The polar coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    DSpanRP operator()(DSpanRP inout_data, ddc::ChunkSpan<CoordRP const, IDomainRP> coordinates)
            const override;
};



/**
 * @brief A class which stores information necessary to create a pointer to an instance of the SplineInterpolatorRP class.
 *
 * This class allows an instance of the SplineInterpolatorRP class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolatorRP to be freed when the object is not in use.
 * These objects are: m_coefs.
 */
class PreallocatableSplineInterpolatorRP : public IPreallocatableInterpolatorRP
{
    SplineRPBuilder const& m_builder;

    SplineRPEvaluator const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolatorRP objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolatorRP(
            SplineRPBuilder const& builder,
            SplineRPEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolatorRP() override = default;

    /**
     * Create a pointer to an instance of the SplineInterpolatorRP class.
     *
     * @return A pointer to an instance of the SplineInterpolatorRP class.
     */
    std::unique_ptr<IInterpolatorRP> preallocate() const override
    {
        return std::make_unique<SplineInterpolatorRP>(m_builder, m_evaluator);
    }
};
