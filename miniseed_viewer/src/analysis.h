// Operations on pairs of traces (difference, cross-correlation). No graphics here.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "data.h"
#include "dsp.h"

namespace msv {

// One gap-free run of samples as 64-bit numbers.
struct Piece {
    double t0 = 0;  // time of the first sample, seconds since 1970 (UTC)
    double dt = 1;  // seconds per sample
    std::vector<double> v;
    double t_end() const { return t0 + dt * (double)(v.empty() ? 0 : v.size() - 1); }
};
using Pieces = std::vector<Piece>;

// Takes the samples of tr that lie between t1 and t2 (after the shift), multiplies them by scale,
// flips the sign if flip is set, and moves them later in time by shift seconds
// (a negative shift moves them earlier). t1 and t2 are times on the shifted axis.
Pieces prepare_pieces(const Trace& tr, double scale, bool flip, double shift, double t1, double t2);

// If the two sets of pieces have different sample rates, brings the faster one down to the slower one's rate:
// zero-phase low-pass filter (4th-order Butterworth at 40% of the slower rate, so nothing folds back),
// then samples at the slower rate's times (linear interpolation; exact when the rates are in a whole ratio).
// Sets note to say what was done ("" if nothing). Returns false (with err) if either has no data.
bool match_rates(Pieces& a, Pieces& b, std::string& note, std::string& err);

// a - b, on the sample times of a. Where b's samples fall between a's, b is interpolated linearly.
// Only times where both traces have data are kept. Returns false and sets err if that is impossible
// (different sample rates, or no time in common).
bool diff_pieces(const Pieces& a, const Pieces& b, std::vector<Segment>& out, std::string& err);

// Same, but with balance set the output is sqrt(R)*a - b/sqrt(R), where R is the factor that makes the standard
// deviation of the output smallest (the output's own mean does not count). That R is std(b)/std(a), both taken
// over the stretch the traces share. R is returned in *R_out.
// If fa and fb are given (filtered copies of a and b, same times), R is found from those instead, but the
// output is still made from the unfiltered a and b.
bool diff_pieces(const Pieces& a, const Pieces& b, bool balance, double* R_out, std::vector<Segment>& out,
                 std::string& err, const Pieces* fa = nullptr, const Pieces* fb = nullptr);

struct XcorrResult {
    std::vector<double> lag;  // seconds
    std::vector<double> c;    // the cross-correlation at each lag
    double peak_lag = 0, peak_value = 0;  // the largest |c|
    int64_t n_samples = 0;                // length of the stretch that was correlated
    double dt = 1;
};

// c(lag) = sum over t of x1(t) * x2(t + lag), after removing each trace's mean. A peak at a positive lag means
// that trace 2 is a delayed copy of trace 1. Time with no data in a trace counts as zero.
// With normalize, c is divided by sqrt(sum(x1^2) * sum(x2^2)), so identical traces give 1 at lag 0.
// Only lags between lag_min and lag_max (seconds) are returned.
bool cross_correlate(const Pieces& a, const Pieces& b, double lag_min, double lag_max, bool normalize,
                     XcorrResult& out, std::string& err);

// Removes each piece's mean, then filters it (at the piece's own sample rate). The caller should have checked the
// filter with check_filter() for that rate.
void filter_pieces(Pieces& p, const FilterSpec& f);

// The samples of p that lie in [t1, t2).
Pieces cut_pieces(const Pieces& p, double t1, double t2);

// Cross-correlation coefficients in a sequence of windows ("correlogram"). Column k holds the normalized
// cross-correlation (as cross_correlate with normalize) of the two traces within the window that starts at
// t_first + k*hop; rows go from the smallest to the largest lag. NaN where a window has no usable data.
struct CorrelogramResult {
    int ncols = 0, nlags = 0;
    std::vector<float> c;  // nlags * ncols, row r = lag lag0 + r*dt
    double t_first = 0;    // start of the first window (seconds since 1970)
    double hop = 0;        // seconds between window starts
    double win = 0;        // window length, seconds
    double lag0 = 0, dt = 1;
};
// The windows cover the time both traces share. overlap_pct is 0..95.
bool correlogram(const Pieces& a, const Pieces& b, double win_s, int overlap_pct, double lag_min, double lag_max,
                 CorrelogramResult& out, std::string& err);

}  // namespace msv
