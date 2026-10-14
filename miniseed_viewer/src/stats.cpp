#include "stats.h"

#include <algorithm>
#include <cmath>

namespace msv {

namespace {

struct Sums {
    double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    void add(double y) {
        double y2 = y * y;
        s1 += y;
        s2 += y2;
        s3 += y2 * y;
        s4 += y2 * y2;
    }
    void remove(double y) {
        double y2 = y * y;
        s1 -= y;
        s2 -= y2;
        s3 -= y2 * y;
        s4 -= y2 * y2;
    }
};

template <class T>
MovingStats run(const T* x, int64_t n, double t0, double dt, double t1, double t2, double win_s, int max_points) {
    MovingStats out;
    const int64_t W = std::max<int64_t>(8, (int64_t)std::llround(win_s / dt));
    if (W > n) return out;
    const int64_t half = W / 2;
    // Centres whose whole window is inside the data and whose time is between t1 and t2.
    int64_t c0 = (int64_t)std::ceil((t1 - t0) / dt - 1e-9), c1 = (int64_t)std::floor((t2 - t0) / dt + 1e-9);
    c0 = std::max(c0, half);
    c1 = std::min(c1, n - W + half);
    if (c1 < c0) return out;
    const int64_t step = std::max<int64_t>(1, (c1 - c0 + 1) / std::max(1, max_points));

    // Work with (x - m) / s, so that the sums stay well-behaved for data with a large offset.
    const int64_t a = c0 - half, b = c1 - half + W;  // all samples that any window touches
    double m = 0;
    for (int64_t i = a; i < b; i++) m += x[i];
    m /= (double)(b - a);
    double v = 0;
    for (int64_t i = a; i < b; i++) v += (x[i] - m) * (x[i] - m);
    double s = std::sqrt(v / (double)(b - a));
    if (!(s > 0)) s = 1;
    auto y = [&](int64_t i) { return ((double)x[i] - m) / s; };

    Sums S;
    int64_t ws = -1;  // start of the window the sums currently hold
    int since_refresh = 0;
    for (int64_t c = c0; c <= c1; c += step) {
        const int64_t start = c - half;
        if (ws < 0 || start - ws >= W || ++since_refresh >= 256) {  // (re)compute from scratch
            S = Sums();
            for (int64_t i = start; i < start + W; i++) S.add(y(i));
            since_refresh = 0;
        } else {  // slide the window
            for (int64_t i = ws; i < start; i++) S.remove(y(i));
            for (int64_t i = ws + W; i < start + W; i++) S.add(y(i));
        }
        ws = start;
        const double Wd = (double)W;
        const double mu = S.s1 / Wd, e2 = S.s2 / Wd, e3 = S.s3 / Wd, e4 = S.s4 / Wd;
        const double var = e2 - mu * mu;
        const double mu3 = e3 - 3 * mu * e2 + 2 * mu * mu * mu;
        const double mu4 = e4 - 4 * mu * e3 + 6 * mu * mu * e2 - 3 * mu * mu * mu * mu;
        const double mean_sq = s * s * e2 + 2 * s * m * mu + m * m;  // mean of x^2
        out.t.push_back(t0 + dt * (double)c);
        out.rms.push_back(std::sqrt(std::max(0.0, mean_sq)));
        if (var > 1e-14) {
            out.skew.push_back(mu3 / std::pow(var, 1.5));
            out.kurt.push_back(mu4 / (var * var));
        } else {
            out.skew.push_back(std::nan(""));
            out.kurt.push_back(std::nan(""));
        }
    }
    return out;
}

}  // namespace

MovingStats moving_stats(const Segment& seg, double t1, double t2, double win_s, int max_points) {
    if (seg.kind == Segment::Kind::Int32)
        return run(seg.i32.data(), seg.size(), seg.t0, seg.dt, t1, t2, win_s, max_points);
    return run(seg.f32.data(), seg.size(), seg.t0, seg.dt, t1, t2, win_s, max_points);
}

}  // namespace msv
