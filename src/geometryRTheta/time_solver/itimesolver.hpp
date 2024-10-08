// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief Base class for the time solvers.
 */
class ITimeSolverRTheta
{
public:
    virtual ~ITimeSolverRTheta() = default;

    /**
     * @brief Solves on @f$ T = dt*N @f$ the equations system.
     *
     * @param[in, out] allfdistribu
     *      On input: the initial condition.
     *      On output: the solution at @f$ dt *N@f$.
     * @param[in] dt
     *      The time step.
     * @param[in] steps
     *      The number @f$ N@f$ of time interations.
     *
     * @return A Field toward allfdistribu.
     */
    virtual host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            double const dt,
            int const steps = 1) const = 0;


protected:
    /**
     * @brief Displays the time difference between two given times and a title.
     *
     * Useful to display the duration of a simulation.
     *
     * @param[in] title
     *      A string printed in front of the computed duration.
     * @param[in] start_time
     *      The start time point.
     * @param[in] end_time
     *      The end time point.
     */
    void display_time_difference(
            std::string const& title,
            std::chrono::time_point<std::chrono::system_clock> const& start_time,
            std::chrono::time_point<std::chrono::system_clock> const& end_time) const
    {
        double const time_h
                = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();
        double const time_min
                = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
        double const time_s
                = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        double const time_ms
                = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time)
                          .count();
        std::cout << title << time_h << "h " << time_min - time_h * 60 << "min "
                  << time_s - time_min * 60 << "s " << time_ms - time_s * 1000 << "ms "
                  << std::endl;
    }
};
