#include "analysis.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>

#include "dsp.h"
#include "parallel.h"
#include "pocketfft_hdronly.h"

namespace msv {

namespace {

using cd = std::complex<double>;

// Times are seconds since 1970 as doubles, good to about 0.2 microseconds. Samples closer than this fraction
// of a sample interval to the grid are taken to be on it, so that shifted traces still line up exactly.
const double kSnap = 1e-3;

bool same_rate(const Pieces& a, const Pieces& b) {
    return std::fabs(a[0].dt - b[0].dt) <= 1e-6 * a[0].dt;
}

// Value of piece p at time t, which must lie within [p.t0, p.t_end()] (give or take a hair).
double value_at(const Piece& p, double t) {
    double u = (t - p.t0) / p.dt;
    double r = std::round(u);
    int64_t n = (int64_t)p.v.size();
    if (std::fabs(u - r) < kSnap) return p.v[(size_t)std::max<int64_t>(0, std::min<int64_t>(n - 1, (int64_t)r))];
    int64_t j = (int64_t)std::floor(u);
    j = std::max<int64_t>(0, std::min<int64_t>(n - 2, j));
    double w = u - (double)j;
    return p.v[(size_t)j] * (1 - w) + p.v[(size_t)j + 1] * w;
}

}  // namespace

bool match_rates(Pieces& a, Pieces& b, std::string& note, std::string& err) {
    note.clear();
    if (a.empty() || b.empty()) {
        err = "One of the traces has no data in the chosen time range.";
        return false;
    }
    if (same_rate(a, b)) return true;
    const bool a_fast = a[0].dt < b[0].dt;
    Pieces& fast = a_fast ? a : b;
    const double slow_dt = a_fast ? b[0].dt : a[0].dt;
    const double fast_rate = 1.0 / fast[0].dt, slow_rate = 1.0 / slow_dt;
    FilterSpec f;
    f.type = FilterType::LowPass;
    f.order = 4;
    f.f1 = 0.4 * slow_rate;
    f.zero_phase = true;
    const bool can_filter = check_filter(f, fast_rate).empty();
    Pieces out;
    for (Piece& p : fast) {
        if (can_filter) apply_filter(f, fast_rate, p.v);
        Piece q;
        q.t0 = p.t0;
        q.dt = slow_dt;
        const double span = p.t_end() - p.t0;
        const int64_t n = (int64_t)std::floor(span / slow_dt + kSnap) + 1;
        q.v.resize((size_t)n);
        for (int64_t k = 0; k < n; k++) q.v[(size_t)k] = value_at(p, p.t0 + slow_dt * (double)k);
        out.push_back(std::move(q));
    }
    fast = std::move(out);
    char buf[160];
    std::snprintf(buf, sizeof buf, "resampled the %.6g/s trace to %.6g/s (low-pass at %.6g Hz)", fast_rate, slow_rate,
                  f.f1);
    note = buf;
    return true;
}

Pieces prepare_pieces(const Trace& tr, double scale, bool flip, double shift, double t1, double t2) {
    Pieces out;
    for (const Segment& raw : extract(tr, t1 - shift, t2 - shift)) {
        Piece p;
        p.t0 = raw.t0 + shift;
        p.dt = raw.dt;
        p.v = raw.doubles();
        const double k = flip ? -scale : scale;
        if (k != 1.0)
            for (double& d : p.v) d *= k;
        out.push_back(std::move(p));
    }
    return out;
}

bool diff_pieces(const Pieces& a, const Pieces& b, std::vector<Segment>& out, std::string& err) {
    return diff_pieces(a, b, false, nullptr, out, err);
}

namespace {

// Pairs of values at the same times, one run per gap-free stretch the two traces share.
struct DiffRun {
    double t0, dt;
    std::vector<double> va, vb;
};

std::vector<DiffRun> build_runs(const Pieces& a, const Pieces& b) {
    std::vector<DiffRun> runs;
    const double eps = kSnap * a[0].dt;
    for (const Piece& pa : a) {
        const int64_t na = (int64_t)pa.v.size();
        for (const Piece& pb : b) {
            double lo = std::max(pa.t0, pb.t0), hi = std::min(pa.t_end(), pb.t_end());
            if (hi < lo - eps) continue;
            int64_t i0 = (int64_t)std::ceil((pb.t0 - eps - pa.t0) / pa.dt);
            int64_t i1 = (int64_t)std::floor((pb.t_end() + eps - pa.t0) / pa.dt);
            i0 = std::max<int64_t>(i0, 0);
            i1 = std::min<int64_t>(i1, na - 1);
            if (i1 < i0) continue;
            DiffRun r;
            r.t0 = pa.t0 + pa.dt * (double)i0;
            r.dt = pa.dt;
            for (int64_t i = i0; i <= i1; i++) {
                r.va.push_back(pa.v[(size_t)i]);
                r.vb.push_back(value_at(pb, pa.t0 + pa.dt * (double)i));
            }
            runs.push_back(std::move(r));
        }
    }
    return runs;
}

}  // namespace

bool diff_pieces(const Pieces& a, const Pieces& b, bool balance, double* R_out, std::vector<Segment>& out,
                 std::string& err, const Pieces* fa, const Pieces* fb) {
    out.clear();
    if (a.empty() || b.empty()) {
        err = "One of the traces has no data in the chosen time range.";
        return false;
    }
    if (!same_rate(a, b)) {
        err = "The traces have different sample rates (" + std::to_string(1.0 / a[0].dt) + " and " +
              std::to_string(1.0 / b[0].dt) + " per second).";
        return false;
    }
    std::vector<DiffRun> runs = build_runs(a, b);
    if (runs.empty()) {
        err = "The two traces have no time in common (after the time shifts).";
        return false;
    }
    double s = 1.0;  // sqrt(R)
    if (balance) {
        // The statistics come from the filtered copies when there are any (same times, so same runs).
        std::vector<DiffRun> fruns;
        if (fa && fb) fruns = build_runs(*fa, *fb);
        const std::vector<DiffRun>& st = fruns.empty() ? runs : fruns;
        double ma = 0, mb = 0;
        int64_t n = 0;
        for (const DiffRun& r : st)
            for (size_t i = 0; i < r.va.size(); i++) {
                ma += r.va[i];
                mb += r.vb[i];
                n++;
            }
        ma /= (double)n;
        mb /= (double)n;
        double saa = 0, sbb = 0;
        for (const DiffRun& r : st)
            for (size_t i = 0; i < r.va.size(); i++) {
                saa += (r.va[i] - ma) * (r.va[i] - ma);
                sbb += (r.vb[i] - mb) * (r.vb[i] - mb);
            }
        if (!(saa > 0) || !(sbb > 0)) {
            err = "One of the traces is constant over the shared time, so no scaling factor can be found.";
            return false;
        }
        const double R = std::sqrt(sbb / saa);
        s = std::sqrt(R);
        if (R_out) *R_out = R;
    }
    for (const DiffRun& r : runs) {
        std::vector<double> d(r.va.size());
        for (size_t i = 0; i < d.size(); i++) d[i] = s * r.va[i] - r.vb[i] / s;
        Segment seg;
        seg.t0 = r.t0;
        seg.dt = r.dt;
        seg.set(d);
        out.push_back(std::move(seg));
    }
    return true;
}

void filter_pieces(Pieces& p, const FilterSpec& f) {
    for (Piece& q : p) {
        if (q.v.empty()) continue;
        double m = 0;
        for (double d : q.v) m += d;
        m /= (double)q.v.size();
        for (double& d : q.v) d -= m;
        apply_filter(f, 1.0 / q.dt, q.v);
    }
}

Pieces cut_pieces(const Pieces& p, double t1, double t2) {
    Pieces out;
    for (const Piece& s : p) {
        const double eps = kSnap * s.dt;
        int64_t i0 = (int64_t)std::ceil((t1 - s.t0) / s.dt - kSnap);
        int64_t i1 = (int64_t)std::ceil((t2 - s.t0) / s.dt - kSnap);  // exclusive
        i0 = std::max<int64_t>(i0, 0);
        i1 = std::min<int64_t>(i1, (int64_t)s.v.size());
        (void)eps;
        if (i1 <= i0) continue;
        Piece q;
        q.t0 = s.t0 + s.dt * (double)i0;
        q.dt = s.dt;
        q.v.assign(s.v.begin() + i0, s.v.begin() + i1);
        out.push_back(std::move(q));
    }
    return out;
}

bool cross_correlate(const Pieces& a, const Pieces& b, double lag_min, double lag_max, bool normalize,
                     XcorrResult& R, std::string& err) {
    R = XcorrResult();
    if (a.empty() || b.empty()) {
        err = "One of the traces has no data in the chosen time range.";
        return false;
    }
    if (!same_rate(a, b)) {
        err = "The traces have different sample rates (" + std::to_string(1.0 / a[0].dt) + " and " +
              std::to_string(1.0 / b[0].dt) + " per second).";
        return false;
    }
    if (!(lag_max > lag_min)) {
        err = "The largest lag must be greater than the smallest lag.";
        return false;
    }
    const double dt = a[0].dt;
    const double T0 = a[0].t0;  // trace 1's first sample defines the grid
    // Grid indices (relative to T0) covered by either trace.
    int64_t g_lo = 0, g_hi = (int64_t)std::floor((a.back().t_end() - T0) / dt + 0.5);
    g_lo = std::min<int64_t>(g_lo, (int64_t)std::ceil((b[0].t0 - T0) / dt - kSnap));
    g_hi = std::max<int64_t>(g_hi, (int64_t)std::floor((b.back().t_end() - T0) / dt + kSnap));
    const int64_t N = g_hi - g_lo + 1;
    const int64_t kmin = (int64_t)std::ceil(lag_min / dt - 1e-9), kmax = (int64_t)std::floor(lag_max / dt + 1e-9);
    if (kmax < kmin) {
        err = "The lag range is shorter than one sample.";
        return false;
    }
    if (N > ((int64_t)1 << 26)) {
        err = "Too many samples (" + std::to_string((long long)N) + ") to correlate in one go. Select a shorter time window on Main.";
        return false;
    }
    const int64_t L = std::max(std::llabs(kmin), std::llabs(kmax));
    if (L >= N) {
        err = "The lags reach beyond the length of the data (" + std::to_string((double)N * dt) + " s).";
        return false;
    }
    R.dt = dt;
    R.n_samples = N;

    // Put both traces on the grid (zero where there is no data), with their means removed.
    std::vector<double> x((size_t)N, 0.0), y((size_t)N, 0.0);
    std::vector<char> vx((size_t)N, 0), vy((size_t)N, 0);
    for (const Piece& p : a)
        for (size_t i = 0; i < p.v.size(); i++) {
            int64_t g = (int64_t)std::llround((p.t0 + dt * (double)i - T0) / dt) - g_lo;
            if (g >= 0 && g < N) {
                x[(size_t)g] = p.v[i];
                vx[(size_t)g] = 1;
            }
        }
    const double eps = kSnap * dt;
    for (const Piece& p : b) {
        int64_t i0 = (int64_t)std::ceil((p.t0 - eps - T0) / dt) - g_lo;
        int64_t i1 = (int64_t)std::floor((p.t_end() + eps - T0) / dt) - g_lo;
        for (int64_t g = std::max<int64_t>(i0, 0); g <= std::min<int64_t>(i1, N - 1); g++) {
            y[(size_t)g] = value_at(p, T0 + dt * (double)(g + g_lo));
            vy[(size_t)g] = 1;
        }
    }
    double e1 = 0, e2 = 0;
    for (auto* pr : {&x, &y}) {
        const std::vector<char>& v = pr == &x ? vx : vy;
        double sum = 0;
        int64_t cnt = 0;
        for (int64_t i = 0; i < N; i++)
            if (v[(size_t)i]) {
                sum += (*pr)[(size_t)i];
                cnt++;
            }
        if (!cnt) {
            err = "One of the traces has no data in the chosen time range.";
            return false;
        }
        const double mean = sum / (double)cnt;
        double e = 0;
        for (int64_t i = 0; i < N; i++) {
            if (v[(size_t)i]) {
                (*pr)[(size_t)i] -= mean;
                e += (*pr)[(size_t)i] * (*pr)[(size_t)i];
            }
        }
        (pr == &x ? e1 : e2) = e;
    }
    if (normalize && (e1 <= 0 || e2 <= 0)) {
        err = "One of the traces is constant, so it cannot be normalized.";
        return false;
    }

    // Circular correlation by FFT; padding by L keeps wrapped-around terms out of the lags we want.
    const size_t M = pocketfft::detail::util::good_size_real((size_t)(N + L));
    std::vector<double> bx(M, 0.0), by(M, 0.0);
    std::copy(x.begin(), x.end(), bx.begin());
    std::copy(y.begin(), y.end(), by.begin());
    x.clear();
    y.clear();
    std::vector<cd> X(M / 2 + 1), Y(M / 2 + 1);
    const pocketfft::shape_t shape{M};
    const pocketfft::stride_t sd{(ptrdiff_t)sizeof(double)}, sc{(ptrdiff_t)sizeof(cd)};
    pocketfft::r2c(shape, sd, sc, 0, true, bx.data(), X.data(), 1.0);
    pocketfft::r2c(shape, sd, sc, 0, true, by.data(), Y.data(), 1.0);
    for (size_t k = 0; k < X.size(); k++) X[k] = std::conj(X[k]) * Y[k];
    pocketfft::c2r(shape, sc, sd, 0, false, X.data(), bx.data(), 1.0 / (double)M);

    const double norm = normalize ? 1.0 / std::sqrt(e1 * e2) : 1.0;
    R.lag.reserve((size_t)(kmax - kmin + 1));
    R.c.reserve((size_t)(kmax - kmin + 1));
    double best = -1;
    for (int64_t k = kmin; k <= kmax; k++) {
        double c = bx[(size_t)(((k % (int64_t)M) + (int64_t)M) % (int64_t)M)] * norm;
        R.lag.push_back((double)k * dt);
        R.c.push_back(c);
        if (std::fabs(c) > best) {
            best = std::fabs(c);
            R.peak_lag = (double)k * dt;
            R.peak_value = c;
        }
    }
    return true;
}

bool correlogram(const Pieces& a, const Pieces& b, double win_s, int overlap_pct, double lag_min, double lag_max,
                 CorrelogramResult& R, std::string& err) {
    R = CorrelogramResult();
    if (a.empty() || b.empty()) {
        err = "One of the traces has no data in the chosen time range.";
        return false;
    }
    if (!same_rate(a, b)) {
        err = "The traces have different sample rates.";
        return false;
    }
    if (!(lag_max > lag_min)) {
        err = "The largest lag must be greater than the smallest lag.";
        return false;
    }
    const double dt = a[0].dt;
    const int64_t kmin = (int64_t)std::ceil(lag_min / dt - 1e-9), kmax = (int64_t)std::floor(lag_max / dt + 1e-9);
    if (kmax < kmin) {
        err = "The lag range is shorter than one sample.";
        return false;
    }
    const int64_t L = std::max(std::llabs(kmin), std::llabs(kmax));
    if (!(win_s > 0) || win_s / dt < 2.0 * (double)L + 8) {
        err = "The window must be longer than twice the largest lag (and than a few samples). It is best to make it "
              "much longer, because the overlap between the traces shrinks as the lag grows.";
        return false;
    }
    // the data cover each sample's whole interval, so the span runs one sample past the last time stamp
    const double t_start = std::max(a.front().t0, b.front().t0),
                 t_end = std::min(a.back().t_end(), b.back().t_end()) + dt;
    if (!(t_end - t_start >= win_s - 1e-9)) {
        err = "The window is longer than the time the two traces share (" + std::to_string(t_end - t_start) + " s).";
        return false;
    }
    overlap_pct = std::max(0, std::min(95, overlap_pct));
    const double hop = win_s * (100 - overlap_pct) / 100.0;
    const int64_t ncols = (int64_t)std::floor((t_end - t_start - win_s) / hop + 1e-9) + 1;
    if (ncols > 6000) {
        err = std::to_string((long long)ncols) + " windows would be needed; use a longer window, less overlap or a shorter interval (at most 6000).";
        return false;
    }
    R.ncols = (int)ncols;
    R.nlags = (int)(kmax - kmin + 1);
    R.t_first = t_start;
    R.hop = hop;
    R.win = win_s;
    R.lag0 = (double)kmin * dt;
    R.dt = dt;
    R.c.assign((size_t)R.ncols * R.nlags, std::nanf(""));
    parallel_for(R.ncols, [&](int k) {
        const double w0 = t_start + hop * k, w1 = w0 + win_s;
        Pieces pa = cut_pieces(a, w0, w1), pb = cut_pieces(b, w0, w1);
        XcorrResult X;
        std::string e;
        // Skip windows that are mostly empty: the coefficients would mean little.
        auto count = [](const Pieces& p) {
            size_t n = 0;
            for (const Piece& q : p) n += q.v.size();
            return n;
        };
        if (pa.empty() || pb.empty() || count(pa) < 0.5 * win_s / dt || count(pb) < 0.5 * win_s / dt) return;
        if (!cross_correlate(pa, pb, R.lag0, (double)kmax * dt, true, X, e)) return;
        if ((int)X.c.size() != R.nlags) return;
        for (int r = 0; r < R.nlags; r++) R.c[(size_t)r * R.ncols + k] = (float)X.c[(size_t)r];
    });
    return true;
}

}  // namespace msv
