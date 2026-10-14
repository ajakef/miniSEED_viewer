// Moving statistics (RMS, skewness, kurtosis) over a time series.
#pragma once
#include <vector>

#include "data.h"

namespace msv {

struct MovingStats {
    std::vector<double> t;     // time (seconds since 1970) of the centre of each window
    std::vector<double> rms;   // root-mean-square amplitude in the window
    std::vector<double> skew;  // skewness: third central moment / variance^1.5 (NaN if flat)
    std::vector<double> kurt;  // kurtosis: fourth central moment / variance^2, 3 for a normal distribution (NaN if flat)
};

// Statistics in windows of win_s seconds centred at about max_points evenly spaced times between
// t1 and t2. Only windows that lie completely inside the segment are used.
MovingStats moving_stats(const Segment& seg, double t1, double t2, double win_s, int max_points);

}  // namespace msv
