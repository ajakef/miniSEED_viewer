// Loading miniSEED files and holding the samples in memory.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace msv {

// Points for drawing a line: x = time (seconds since 1970, UTC), y = amplitude.
struct Series {
    std::vector<double> x, y;
    void clear() { x.clear(); y.clear(); }
};

// A run of evenly spaced samples with no gaps. Samples from a file are kept as
// 32-bit integers (exactly as recorded); results of filtering etc. are kept as
// 32-bit floats. Either way, use at() / doubles() to read them.
struct Segment {
    enum class Kind { Int32, Float32 };

    double t0 = 0;  // time of first sample, seconds since 1970 (UTC)
    double dt = 1;  // seconds per sample
    Kind kind = Kind::Float32;
    std::vector<int32_t> i32;
    std::vector<float> f32;

    int64_t size() const { return (int64_t)(kind == Kind::Int32 ? i32.size() : f32.size()); }
    double at(int64_t i) const { return kind == Kind::Int32 ? (double)i32[i] : (double)f32[i]; }
    double t_end() const { return t0 + dt * (double)(size() ? size() - 1 : 0); }
    double time_at(int64_t i) const { return t0 + dt * (double)i; }

    // Samples [i0, i1) as 64-bit numbers.
    std::vector<double> doubles(int64_t i0, int64_t i1) const;
    std::vector<double> doubles() const { return doubles(0, size()); }
    // Replaces the contents with v (stored as 32-bit floats).
    void set(const std::vector<double>& v);
    // Appends samples [skip, size) of other, which must have the same spacing.
    void append(const Segment& other, int64_t skip = 0);

    // Index range [i0, i1) of the samples between times t1 and t2.
    void index_range(double t1, double t2, int64_t& i0, int64_t& i1) const;

    // Summaries of min/max over blocks of samples, so that zoomed-out plots
    // do not have to read every sample. Call build_summary() after the samples are final.
    void build_summary();
    void minmax(int64_t i0, int64_t i1, double& mn, double& mx) const;

    // Fills out with points to draw between times t1 and t2 using about nbins
    // pixel columns. Zoomed far out: a min and a max point per column.
    // Zoomed far in: the samples themselves.
    void envelope(double t1, double t2, int nbins, Series& out) const;

    // Root-mean-square of samples between t1 and t2 (NaN if there are none).
    double rms(double t1, double t2) const;

private:
    struct Level {
        int64_t block;
        std::vector<double> mn, mx;
    };
    std::vector<Level> levels_;
};

// One channel: all the segments found for one network.station.location.channel.
struct Trace {
    std::string id;  // NET.STA.LOC.CHA
    std::string net, sta, loc, cha;
    std::string file;      // first file this channel was found in
    int file_order = 0;    // position of that file on the command line
    std::vector<Segment> segs;
    double t_min = 0, t_max = 0;  // time of first and last sample
    double v_min = 0, v_max = 0;  // smallest and largest sample
    double rate = 0;              // samples per second (first segment)
    bool on = true;               // shown in the plots
};

struct LoadResult {
    std::vector<Trace> traces;
    std::vector<std::string> messages;  // problems, one per line
};

// Reads miniSEED files (version 2 or 3, any encoding libmseed understands) and
// merges them into traces, in order of first appearance.
LoadResult load_files(const std::vector<std::string>& files);

// Appends the contents of more files to an existing set of traces.
void add_files(std::vector<Trace>& traces, const std::vector<std::string>& files,
               std::vector<std::string>& messages);

// Copies the samples of tr between t1 and t2 (one Segment per gap-free run).
// The copies are stored in whatever form the trace uses.
std::vector<Segment> extract(const Trace& tr, double t1, double t2);

// "2026-10-10 03:15:26.510" (UTC).
std::string format_time(double t, int decimals = 3);

}  // namespace msv
