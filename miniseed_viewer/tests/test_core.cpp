// Checks of the non-graphical code. Compared against ObsPy/SciPy by tests/check_core.py.
//   test_core load FILE...       prints a summary of what was read
//   test_core dump OUTDIR        writes synthetic data, filtered data and spectra to OUTDIR
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "../src/analysis.h"
#include "../src/data.h"
#include "../src/dsp.h"
#include "../src/picks.h"
#include "../src/stats.h"

using namespace msv;

static void write_vec(const std::string& path, const std::vector<double>& v) {
    FILE* f = fopen(path.c_str(), "wb");
    fwrite(v.data(), sizeof(double), v.size(), f);
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc >= 2 && !strcmp(argv[1], "analysis")) {  // difference and cross-correlation on synthetic data
        int bad = 0;
        auto check = [&](bool ok, const char* what) {
            printf("%s %s\n", ok ? "ok  " : "FAIL", what);
            if (!ok) bad++;
        };
        const double rate = 100, dt = 1 / rate, t0 = 1.6e9;
        const int n = 6000;
        // noise-like signal; trace 2 is trace 1 delayed by 0.37 s (37 samples), scaled by 2
        std::vector<double> sig(n + 100);
        unsigned seed = 12345;
        for (double& v : sig) {
            seed = seed * 1664525u + 1013904223u;
            v = (double)(seed >> 8) / 16777216.0 - 0.5;
        }
        Trace a, b;
        Segment sa, sb;
        sa.t0 = sb.t0 = t0;
        sa.dt = sb.dt = dt;
        sa.kind = sb.kind = Segment::Kind::Int32;
        for (int i = 0; i < n; i++) {
            sa.i32.push_back((int)std::lround(sig[i + 50] * 1e5));
            sb.i32.push_back((int)std::lround(2 * sig[i + 50 - 37] * 1e5));
        }
        a.segs.push_back(sa);
        b.segs.push_back(sb);
        Pieces pa = prepare_pieces(a, 1.0, false, 0, -1e300, 1e300);
        Pieces pb = prepare_pieces(b, 1.0, false, 0, -1e300, 1e300);
        XcorrResult R;
        std::string err;
        check(cross_correlate(pa, pb, -1, 1, true, R, err), "xcorr runs");
        check(R.lag.size() == 201, "xcorr returns 201 lags for -1..1 s");
        check(std::fabs(R.peak_lag - 0.37) < 1e-9, "peak at +0.37 s (trace 2 later)");
        check(R.peak_value > 0.95, "peak value close to 1");
        // shifting trace 2 earlier by 0.37 s should put the peak at 0
        Pieces pb2 = prepare_pieces(b, 1.0, false, -0.37, -1e300, 1e300);
        cross_correlate(pa, pb2, -1, 1, true, R, err);
        check(std::fabs(R.peak_lag) < 1e-9, "after shifting trace 2 by -0.37 s the peak is at 0");
        // unnormalized, versus brute force at lag +37
        cross_correlate(pa, pb, 0.365, 0.375, false, R, err);
        double mx = 0, my = 0;
        for (int i = 0; i < n; i++) {
            mx += sa.i32[i];
            my += sb.i32[i];
        }
        mx /= n;
        my /= n;
        double brute = 0;
        for (int i = 0; i + 37 < n; i++) brute += (sa.i32[i] - mx) * (sb.i32[i + 37] - my);
        check(R.c.size() == 1 && std::fabs(R.c[0] - brute) < 1e-6 * std::fabs(brute), "unnormalized value matches a direct sum");
        // negative lags
        cross_correlate(pb, pa, -1, 1, true, R, err);
        check(std::fabs(R.peak_lag + 0.37) < 1e-9, "swapping the traces mirrors the lag");
        // difference: a - (b/2 shifted earlier by 0.37 s) is small (only rounding)
        Trace b2;
        Segment sh = sb;
        for (auto& v : sh.i32) v /= 2;
        b2.segs.push_back(sh);
        std::vector<Segment> out;
        Pieces pbs = prepare_pieces(b2, 1.0, false, -0.37, -1e300, 1e300);
        check(diff_pieces(pa, pbs, out, err), "diff runs");
        check(out.size() == 1 && out[0].size() == n - 37, "diff keeps only the overlapping n-37 samples");
        double worst = 0;
        for (int64_t i = 0; i < out[0].size(); i++) worst = std::max(worst, std::fabs(out[0].at(i)));
        check(worst <= 1.0, "diff is zero apart from integer rounding");
        // half-sample shift: result interpolated, grid stays trace 1's
        Pieces ph = prepare_pieces(a, 1.0, false, 0.005, -1e300, 1e300);
        check(diff_pieces(pa, ph, out, err) && out[0].size() == n - 1 && std::fabs(out[0].t0 - (t0 + dt)) < 1e-9,
              "half-sample shift: n-1 samples on trace 1's grid");
        // different rates are refused; no overlap is reported
        Pieces pr = pb;
        pr[0].dt = 0.005;
        check(!diff_pieces(pa, pr, out, err), "different rates refused");
        Pieces far = prepare_pieces(b, 1.0, false, 1000, -1e300, 1e300);
        check(!diff_pieces(pa, far, out, err), "no overlap reported");
        // window: only the middle 10 s
        Pieces pw = prepare_pieces(a, 1.0, false, 0, t0 + 10, t0 + 20);
        check(pw.size() == 1 && pw[0].v.size() == 1001, "window 10..20 s holds 1001 samples");
        // balanced difference: trace 2 = 3 * trace 1 (plus an offset) must give R = 3 and a flat output
        {
            Trace c3;
            Segment s3 = sa;
            for (auto& v : s3.i32) v = 3 * v + 1000;
            c3.segs.push_back(s3);
            Pieces p3 = prepare_pieces(c3, 1.0, false, 0, -1e300, 1e300);
            double R = 0;
            check(diff_pieces(pa, p3, true, &R, out, err), "balanced diff runs");
            check(std::fabs(R - 3) < 1e-9, "R = std(tr2)/std(tr1) = 3");
            double m = 0, w = 0;
            for (int64_t i = 0; i < out[0].size(); i++) m += out[0].at(i);
            m /= (double)out[0].size();
            for (int64_t i = 0; i < out[0].size(); i++) w = std::max(w, std::fabs(out[0].at(i) - m));
            check(w < 1e-3 * 1e5, "balanced output is flat (apart from its mean)");
            // R really is the minimum: the std of the output grows for R*1.1 and R/1.1
            auto sd = [&](double Rr) {
                double s = std::sqrt(Rr), m1 = 0, v2 = 0;
                int n1 = (int)pa[0].v.size();
                std::vector<double> d(n1);
                for (int i = 0; i < n1; i++) { d[i] = s * pa[0].v[i] - p3[0].v[i] / s; m1 += d[i]; }
                m1 /= n1;
                for (double x : d) v2 += (x - m1) * (x - m1);
                return v2;
            };
            check(sd(3) < sd(3.3) && sd(3) < sd(3 / 1.1), "R is the minimum of the variance");
        }
        // R from filtered copies: trace 2 = 3 * trace 1 + a huge slow drift; output stays unfiltered
        {
            Trace d3;
            Segment s3 = sa;
            for (int i = 0; i < (int)s3.i32.size(); i++) s3.i32[i] = 3 * s3.i32[i] + (int)(5e6 * std::sin(2 * M_PI * 0.01 * i * dt));
            d3.segs.push_back(s3);
            Pieces p3 = prepare_pieces(d3, 1.0, false, 0, -1e300, 1e300);
            double R0 = 0, R1 = 0;
            diff_pieces(pa, p3, true, &R0, out, err);
            FilterSpec hp;
            hp.type = FilterType::HighPass;
            hp.f1 = 5;
            hp.order = 4;
            Pieces fa = pa, fb = p3;
            for (Pieces* ps : {&fa, &fb})
                for (Piece& p : *ps) apply_filter(hp, 1.0 / p.dt, p.v);
            check(diff_pieces(pa, p3, true, &R1, out, err, &fa, &fb), "balanced diff with filtered R runs");
            check(std::fabs(R0 - 3) > 1, "without a filter the drift spoils R");
            check(std::fabs(R1 - 3) < 0.05, "with a high-pass R is about 3");
            double s = std::sqrt(R1);
            check(std::fabs(out[0].at(1000) - (s * pa[0].v[1000] - p3[0].v[1000] / s)) < 1e-2,
                  "output is made from the unfiltered traces");
        }
        // correlogram: trace 2 = trace 1 delayed by 0.37 s, 60 s of noise at 100/s, 20 s windows, 50% overlap
        {
            Piece x1, x2;
            x1.t0 = x2.t0 = t0;
            x1.dt = x2.dt = dt;
            unsigned sd = 777;
            std::vector<double> nz(6100);
            for (double& v : nz) {
                sd = sd * 1664525u + 1013904223u;
                v = (double)(sd >> 8) / 16777216.0 - 0.5;
            }
            for (int i = 0; i < 6000; i++) {
                x1.v.push_back(nz[i + 100]);
                x2.v.push_back(nz[i + 100 - 37]);
            }
            CorrelogramResult C;
            check(correlogram({x1}, {x2}, 20, 50, -1, 1, C, err), "correlogram runs");
            check(C.ncols == 5 && C.nlags == 201, "5 windows (0,10,20,30,40 s) and 201 lags");
            bool ok = true;
            for (int k = 0; k < C.ncols; k++) {
                int best = 0;
                for (int r = 0; r < C.nlags; r++)
                    if (C.c[(size_t)r * C.ncols + k] > C.c[(size_t)best * C.ncols + k]) best = r;
                if (std::fabs(C.lag0 + best * C.dt - 0.37) > 1e-9 || C.c[(size_t)best * C.ncols + k] < 0.9) ok = false;
            }
            check(ok, "every window peaks at +0.37 s with coefficient > 0.9");
            check(!correlogram({x1}, {x2}, 1.5, 50, -1, 1, C, err), "window shorter than twice the lag is refused");
            check(!correlogram({x1}, {x2}, 100, 50, -1, 1, C, err), "window longer than the data is refused");
            Piece g1 = x1;  // a hole in trace 1: windows lying in it come out empty (NaN)
            for (int i = 2000; i < 4000; i++) g1.v[i] = 0;
            Piece ga = x1, gb = x1;
            ga.v.resize(2000);
            gb.t0 = t0 + 40.0;
            gb.v.assign(x1.v.begin() + 4000, x1.v.end());
            Piece x2g = x2;
            CorrelogramResult C2;
            check(correlogram({ga, gb}, {x2g}, 10, 0, -0.5, 0.5, C2, err), "correlogram with a gap runs");
            int nan_cols = 0;
            for (int k = 0; k < C2.ncols; k++)
                if (std::isnan(C2.c[(size_t)(C2.nlags / 2) * C2.ncols + k])) nan_cols++;
            check(nan_cols >= 2, "windows inside the gap are NaN");
        }
        // different rates: 200/s trace with 10 Hz + 90 Hz (the 90 Hz would fold onto 10 Hz at 100/s)
        {
            Trace f, s;
            Segment sf, ss;
            sf.t0 = ss.t0 = t0;
            sf.dt = 0.005;
            ss.dt = 0.01;
            sf.kind = ss.kind = Segment::Kind::Int32;
            for (int i = 0; i < 20000; i++) {
                double t = i * 0.005;
                sf.i32.push_back((int)std::lround(1e4 * (std::sin(2 * M_PI * 10 * t) + std::sin(2 * M_PI * 90 * t))));
            }
            for (int i = 0; i < 10000; i++) ss.i32.push_back((int)std::lround(1e4 * std::sin(2 * M_PI * 10 * i * 0.01)));
            f.segs.push_back(sf);
            s.segs.push_back(ss);
            Pieces pf = prepare_pieces(f, 1.0, false, 0, -1e300, 1e300);
            Pieces ps = prepare_pieces(s, 1.0, false, 0, -1e300, 1e300);
            std::string note;
            check(match_rates(pf, ps, note, err) && !note.empty(), "rates matched, with a note");
            check(std::fabs(pf[0].dt - 0.01) < 1e-12 && pf[0].v.size() == 10000, "fast trace now at 100/s");
            check(diff_pieces(pf, ps, out, err), "diff of mixed rates runs");
            double w = 0;
            for (int64_t i = 500; i < out[0].size() - 500; i++) w = std::max(w, std::fabs(out[0].at(i)));
            check(w < 300, "90 Hz removed, 10 Hz kept (middle of the record)");
            cross_correlate(pf, ps, -0.5, 0.5, true, R, err);
            check(std::fabs(R.peak_lag) < 1e-9 && R.peak_value > 0.99, "mixed-rate xcorr peaks at lag 0");
        }
        printf("%d failure(s)\n", bad);
        return bad ? 1 : 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "load")) {
        std::vector<std::string> files(argv + 2, argv + argc);
        LoadResult r = load_files(files);
        for (auto& m : r.messages) printf("MSG %s\n", m.c_str());
        for (auto& t : r.traces) {
            printf("TRACE %s nseg=%zu rate=%.6f t_min=%.6f t_max=%.6f v=[%.0f,%.0f]\n", t.id.c_str(), t.segs.size(), t.rate,
                   t.t_min, t.t_max, t.v_min, t.v_max);
            for (auto& s : t.segs) {
                printf("  SEG t0=%.6f n=%lld first=%.0f,%.0f,%.0f last=%.0f rms=%.6f start=%s\n", s.t0, (long long)s.size(),
                       s.at(0), s.at(1), s.at(2), s.at(s.size() - 1), s.rms(s.t0, s.t_end()), format_time(s.t0).c_str());
            }
        }
        // Check the min/max summary against a plain scan for some ranges.
        for (auto& t : r.traces)
            for (auto& s : t.segs) {
                int bad = 0;
                for (int k = 0; k < 200; k++) {
                    int64_t a = (k * 7919) % (s.size() / 2), b = a + 1 + (k * 104729) % (s.size() - a - 1);
                    double mn, mx, mn2 = 1e300, mx2 = -1e300;
                    s.minmax(a, b, mn, mx);
                    for (int64_t i = a; i < b; i++) {
                        mn2 = std::fmin(mn2, s.at(i));
                        mx2 = std::fmax(mx2, s.at(i));
                    }
                    if (mn != mn2 || mx != mx2) bad++;
                }
                printf("  MINMAX_BAD %d\n", bad);
                Series ser;
                s.envelope(s.t0 + 100, s.t0 + 400, 1000, ser);
                printf("  ENV zoomed-in points=%zu\n", ser.x.size());
                s.envelope(s.t0, s.t_end(), 1500, ser);
                printf("  ENV full points=%zu\n", ser.x.size());
            }
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "time")) {  // how long does loading take?
        std::vector<std::string> files(argv + 2, argv + argc);
        std::vector<Trace> traces;
        std::vector<std::string> msgs;
        auto t0 = std::chrono::steady_clock::now();
        scan_files(traces, files, msgs);
        auto t1 = std::chrono::steady_clock::now();
        std::vector<Trace*> all;
        for (auto& t : traces) all.push_back(&t);
        load_traces(all, msgs);
        auto t2 = std::chrono::steady_clock::now();
        int64_t n = 0;
        for (auto& t : traces)
            for (auto& sg : t.segs) n += sg.size();
        printf("%zu files -> %zu traces, %lld samples: scan %.2f s, load %.2f s\n", files.size(), traces.size(),
               (long long)n, std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count());
        for (auto& m : msgs) printf("MSG %s\n", m.c_str());
        return 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "dump")) {
        std::string out = argv[2];
        const double fs = 200;
        const int n = 20000;
        std::vector<double> x(n);
        unsigned s = 12345;
        for (int i = 0; i < n; i++) {
            s = s * 1664525u + 1013904223u;
            double noise = ((s >> 8) / 16777216.0 - 0.5);
            double t = i / fs;
            x[i] = 1000 + 5 * t + 100 * std::sin(2 * M_PI * 0.3 * t) + 50 * std::sin(2 * M_PI * 5 * t) +
                   20 * std::sin(2 * M_PI * 40 * t) + 10 * noise;
        }
        write_vec(out + "/x.bin", x);
        struct Case { FilterType t; int order; double f1, f2; bool zp; };
        Case cases[] = {{FilterType::LowPass, 4, 10, 0, true},  {FilterType::HighPass, 3, 1, 0, true},
                        {FilterType::BandPass, 4, 1, 10, true}, {FilterType::BandStop, 2, 4, 6, true},
                        {FilterType::BandPass, 3, 2, 20, false}, {FilterType::LowPass, 1, 20, 0, false},
                        {FilterType::BandStop, 3, 4, 40, true}, {FilterType::LowPass, 5, 3, 0, true},
                        {FilterType::BandPass, 5, 0.5, 30, true}, {FilterType::HighPass, 4, 0.1, 0, false}};
        int ci = 0;
        for (auto& c : cases) {
            FilterSpec f;
            f.type = c.t; f.order = c.order; f.f1 = c.f1; f.f2 = c.f2; f.zero_phase = c.zp;
            std::vector<double> y = x;
            apply_filter(f, fs, y);
            write_vec(out + "/y" + std::to_string(ci++) + ".bin", y);
        }
        std::vector<double> d = x;
        detrend(d, Detrend::Linear);
        write_vec(out + "/detrend_linear.bin", d);
        Spectrum p1 = power_spectrum(x.data(), x.size(), fs, false, 0);
        write_vec(out + "/spec1_f.bin", p1.f);
        write_vec(out + "/spec1_p.bin", p1.p);
        Spectrum p2 = power_spectrum(x.data(), x.size(), fs, true, 1024);
        write_vec(out + "/spec2_f.bin", p2.f);
        write_vec(out + "/spec2_p.bin", p2.p);
        {   // moving statistics of the same data, as an integer record (like real data)
            Segment seg;
            seg.kind = Segment::Kind::Int32;
            seg.t0 = 100.0;
            seg.dt = 1.0 / fs;
            for (double v : x) seg.i32.push_back((int32_t)std::lround(v));
            MovingStats ms = moving_stats(seg, 110.0, 190.0, 5.0, 300);
            write_vec(out + "/ms_t.bin", ms.t);
            write_vec(out + "/ms_rms.bin", ms.rms);
            write_vec(out + "/ms_skew.bin", ms.skew);
            write_vec(out + "/ms_kurt.bin", ms.kurt);
        }
        {   // picks file
            std::vector<Pick> picks = {{1791597933.25, ".10024..HDF", "P"}, {1791590400.0, "XX.STA.00.BHZ", "a,b"},
                                       {1791600000.123456, "XX.STA.00.BHZ", "S"}};
            std::string err;
            save_picks(out + "/picks.csv", picks, err);
            std::vector<Pick> back;
            load_picks(out + "/picks.csv", back, err);
            printf("picks round trip: %zu read back, first %s\n", back.size(), format_iso(back[0].t).c_str());
            double tt = 0;
            printf("parse_iso ok=%d t=%.6f\n", (int)parse_iso("2026-10-10T01:25:33.250Z", tt), tt);
        }
        printf("time fmt: %s | %s\n", format_time(1791590126.51).c_str(), format_time(0.9996, 3).c_str());
        return 0;
    }
    fprintf(stderr, "usage: test_core load FILE... | dump OUTDIR\n");
    return 1;
}
