

# File itimesolver.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**time\_solver**](dir_4c2664fc2adc717d896afdb0f76e6fe5.md) **>** [**itimesolver.hpp**](geometryRTheta_2time__solver_2itimesolver_8hpp.md)

[Go to the documentation of this file](geometryRTheta_2time__solver_2itimesolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class ITimeSolverRTheta
{
public:
    virtual ~ITimeSolverRTheta() = default;

    virtual host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> allfdistribu,
            double const dt,
            int const steps = 1) const = 0;


protected:
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
```


