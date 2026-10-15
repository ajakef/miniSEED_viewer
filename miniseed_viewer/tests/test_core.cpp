// Checks of the non-graphical code. Compared against ObsPy/SciPy by tests/check_core.py.
//   test_core load FILE...       prints a summary of what was read
//   test_core dump OUTDIR        writes synthetic data, filtered data and spectra to OUTDIR
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

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
