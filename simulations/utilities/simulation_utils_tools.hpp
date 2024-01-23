
#pragma once

#include <chrono>
#include <fstream>
#include <iostream>



void display_time_difference(
        std::string const& title,
        std::chrono::time_point<std::chrono::system_clock> const& start_time,
        std::chrono::time_point<std::chrono::system_clock> const& end_time)
{
    double const time_h
            = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();
    double const time_min
            = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
    double const time_s
            = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    double const time_ms
            = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << title << time_h << "h " << time_min - time_h * 60 << "min "
              << time_s - time_min * 60 << "s " << time_ms - time_s * 1000 << "ms " << std::endl;
};
