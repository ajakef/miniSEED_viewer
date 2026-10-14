#include "data.h"

#include "parallel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <set>

extern "C" {
#include "libmseed.h"
}

namespace msv {

// ---------------------------------------------------------------------------
// Segment

void Segment::index_range(double t1, double t2, int64_t& i0, int64_t& i1) const {
    if (size() == 0 || t2 < t1) {
        i0 = i1 = 0;
        return;
    }
    double a = std::ceil((t1 - t0) / dt - 1e-6);
    double b = std::floor((t2 - t0) / dt + 1e-6);
    a = std::max(a, 0.0);
    b = std::min(b, (double)(size() - 1));
    if (b < a) {
        i0 = i1 = 0;
        return;
    }
    i0 = (int64_t)a;
    i1 = (int64_t)b + 1;
}

static const int64_t kFirstBlock = 16;
static const int64_t kBlockGrowth = 4;

namespace {
template <class T>
void min_max_raw(const T* p, int64_t i0, int64_t i1, double& mn, double& mx) {
    T lo = p[i0], hi = p[i0];
    for (int64_t i = i0 + 1; i < i1; i++) {
        lo = std::min(lo, p[i]);
        hi = std::max(hi, p[i]);
    }
    mn = std::min(mn, (double)lo);
    mx = std::max(mx, (double)hi);
}
void min_max_raw(const Segment& s, int64_t i0, int64_t i1, double& mn, double& mx) {
    if (i1 <= i0) return;
    if (s.kind == Segment::Kind::Int32)
        min_max_raw(s.i32.data(), i0, i1, mn, mx);
    else
        min_max_raw(s.f32.data(), i0, i1, mn, mx);
}
}  // namespace

void Segment::build_summary() {
    levels_.clear();
    const int64_t n = size();
    int64_t block = kFirstBlock;
    while (block < n) {
        Level L;
        L.block = block;
        int64_t nb = (n + block - 1) / block;
        L.mn.resize(nb);
        L.mx.resize(nb);
        const bool from_raw = levels_.empty();
        auto fill = [&](int64_t b0, int64_t b1) {
            for (int64_t b = b0; b < b1; b++) {
                double lo = std::numeric_limits<double>::infinity(), hi = -lo;
                if (from_raw) {
                    min_max_raw(*this, b * block, std::min(n, (b + 1) * block), lo, hi);
                } else {
                    const Level& P = levels_.back();
                    int64_t s = b * kBlockGrowth, e = std::min<int64_t>((int64_t)P.mn.size(), s + kBlockGrowth);
                    for (int64_t i = s; i < e; i++) {
                        lo = std::min(lo, P.mn[i]);
                        hi = std::max(hi, P.mx[i]);
                    }
                }
                L.mn[b] = lo;
                L.mx[b] = hi;
            }
        };
        if (from_raw && n > 4000000) {  // the first level reads every sample: split it up
            const int chunks = 32;
            parallel_for(chunks, [&](int c) { fill(nb * c / chunks, nb * (c + 1) / chunks); });
        } else {
            fill(0, nb);
        }
        levels_.push_back(std::move(L));
        block *= kBlockGrowth;
    }
}

void Segment::minmax(int64_t i0, int64_t i1, double& mn, double& mx) const {
    mn = std::numeric_limits<double>::infinity();
    mx = -mn;
    i0 = std::max<int64_t>(i0, 0);
    i1 = std::min(i1, size());
    if (i1 <= i0) return;
    // Biggest block size that is at most 1/8 of the range.
    const Level* L = nullptr;
    for (const Level& lv : levels_)
        if (lv.block * 8 <= i1 - i0) L = &lv;
    if (!L) {
        min_max_raw(*this, i0, i1, mn, mx);
        return;
    }
    const int64_t B = L->block;
    int64_t b0 = (i0 + B - 1) / B;  // first whole block
    int64_t b1 = i1 / B;            // one past last whole block
    min_max_raw(*this, i0, std::min(i1, b0 * B), mn, mx);
    for (int64_t b = b0; b < b1; b++) {
        mn = std::min(mn, L->mn[b]);
        mx = std::max(mx, L->mx[b]);
    }
    min_max_raw(*this, std::max(i0, b1 * B), i1, mn, mx);
}

std::vector<double> Segment::doubles(int64_t i0, int64_t i1) const {
    i0 = std::max<int64_t>(i0, 0);
    i1 = std::min(i1, size());
    std::vector<double> v;
    if (i1 <= i0) return v;
    v.resize(i1 - i0);
    if (kind == Kind::Int32)
        for (int64_t i = i0; i < i1; i++) v[i - i0] = i32[i];
    else
        for (int64_t i = i0; i < i1; i++) v[i - i0] = f32[i];
    return v;
}

void Segment::set(const std::vector<double>& v) {
    kind = Kind::Float32;
    i32.clear();
    i32.shrink_to_fit();
    f32.resize(v.size());
    for (size_t i = 0; i < v.size(); i++) f32[i] = (float)v[i];
}

void Segment::append(const Segment& o, int64_t skip) {
    if (kind == o.kind) {
        if (kind == Kind::Int32)
            i32.insert(i32.end(), o.i32.begin() + skip, o.i32.end());
        else
            f32.insert(f32.end(), o.f32.begin() + skip, o.f32.end());
        return;
    }
    // Mixed integer and float data: keep everything as floats.
    if (kind == Kind::Int32) {
        f32.assign(i32.begin(), i32.end());
        i32.clear();
        i32.shrink_to_fit();
        kind = Kind::Float32;
    }
    for (int64_t i = skip; i < o.size(); i++) f32.push_back((float)o.at(i));
}

void Segment::envelope(double t1, double t2, int nbins, Series& out) const {
    out.clear();
    if (size() == 0 || nbins < 1) return;
    // Include one sample beyond each edge so lines reach the plot border.
    int64_t i0 = (int64_t)std::floor((t1 - t0) / dt) - 1;
    int64_t i1 = (int64_t)std::ceil((t2 - t0) / dt) + 2;
    i0 = std::max<int64_t>(i0, 0);
    i1 = std::min(i1, size());
    if (i1 <= i0) return;
    const int64_t cnt = i1 - i0;
    if (cnt <= 2 * (int64_t)nbins) {
        out.x.resize(cnt);
        out.y.resize(cnt);
        for (int64_t i = 0; i < cnt; i++) {
            out.x[i] = time_at(i0 + i);
            out.y[i] = at(i0 + i);
        }
        return;
    }
    out.x.reserve(2 * nbins);
    out.y.reserve(2 * nbins);
    for (int b = 0; b < nbins; b++) {
        int64_t s = i0 + cnt * b / nbins;
        int64_t e = i0 + cnt * (b + 1) / nbins;
        if (e <= s) continue;
        double mn, mx;
        minmax(s, e, mn, mx);
        double tc = t0 + dt * 0.5 * (double)(s + e - 1);
        out.x.push_back(tc);
        out.y.push_back(mn);
        out.x.push_back(tc);
        out.y.push_back(mx);
    }
}

double Segment::rms(double t1, double t2) const {
    int64_t i0, i1;
    index_range(t1, t2, i0, i1);
    if (i1 <= i0) return std::nan("");
    double s = 0;
    for (int64_t i = i0; i < i1; i++) {
        double v = at(i);
        s += v * v;
    }
    return std::sqrt(s / (double)(i1 - i0));
}

// ---------------------------------------------------------------------------
// Loading

namespace {

void parse_sid(const char* sid, Trace& t) {
    char net[64] = "", sta[64] = "", loc[64] = "", cha[64] = "";
    if (ms_sid2nslc((char*)sid, net, sta, loc, cha) == 0) {
        t.net = net;
        t.sta = sta;
        t.loc = loc;
        t.cha = cha;
        t.id = t.net + "." + t.sta + "." + t.loc + "." + t.cha;
    } else {
        t.sta = t.id = sid;
    }
}

bool convert_samples(const MS3TraceSeg* s, Segment& seg) {
    const int64_t n = s->numsamples;
    switch (s->sampletype) {
        case 'i':
            seg.kind = Segment::Kind::Int32;
            seg.i32.assign((const int32_t*)s->datasamples, (const int32_t*)s->datasamples + n);
            return true;
        case 'f':
            seg.kind = Segment::Kind::Float32;
            seg.f32.assign((const float*)s->datasamples, (const float*)s->datasamples + n);
            return true;
        case 'd': {  // rare; stored as 32-bit floats
            const double* p = (const double*)s->datasamples;
            seg.kind = Segment::Kind::Float32;
            seg.f32.resize(n);
            for (int64_t i = 0; i < n; i++) seg.f32[i] = (float)p[i];
            return true;
        }
        default:
            return false;  // text records
    }
}

void finish_trace(Trace& tr) {
    std::sort(tr.segs.begin(), tr.segs.end(),
              [](const Segment& a, const Segment& b) { return a.t0 < b.t0; });
    // First find how big each joined segment will be, so memory is reserved once and exactly.
    std::vector<int64_t> group_size;
    {
        double g_t0 = 0, g_dt = 1;
        for (const Segment& s : tr.segs) {
            if (s.size() == 0) continue;
            if (!group_size.empty() && std::fabs(g_dt - s.dt) < 1e-3 * g_dt) {
                double expect = g_t0 + g_dt * (double)group_size.back();
                if (std::fabs(s.t0 - expect) < 0.5 * g_dt) {
                    group_size.back() += s.size();
                    continue;
                }
                if (s.t0 < expect) {
                    int64_t skip = (int64_t)std::llround((expect - s.t0) / g_dt);
                    if (skip < s.size()) group_size.back() += s.size() - skip;
                    continue;
                }
            }
            group_size.push_back(s.size());
            g_t0 = s.t0;
            g_dt = s.dt;
        }
    }
    // Join segments that touch (within half a sample); trim overlaps.
    std::vector<Segment> merged;
    for (Segment& s : tr.segs) {
        if (s.size() == 0) continue;
        if (!merged.empty()) {
            Segment& p = merged.back();
            bool same_rate = std::fabs(p.dt - s.dt) < 1e-3 * p.dt;
            double expect = p.t_end() + p.dt;
            if (same_rate && std::fabs(s.t0 - expect) < 0.5 * p.dt) {
                p.append(s);
                continue;
            }
            if (same_rate && s.t0 < expect) {  // overlap: keep earlier data
                int64_t skip = (int64_t)std::llround((expect - s.t0) / p.dt);
                if (skip >= s.size()) continue;
                p.append(s, skip);
                continue;
            }
        }
        merged.push_back(std::move(s));
        Segment& b = merged.back();  // room for everything that will be joined on
        const size_t want = (size_t)group_size[merged.size() - 1];
        if (b.kind == Segment::Kind::Int32) b.i32.reserve(want);
        else b.f32.reserve(want);
    }
    tr.segs = std::move(merged);
    tr.t_min = tr.t_max = 0;
    tr.v_min = std::numeric_limits<double>::infinity();
    tr.v_max = -tr.v_min;
    tr.nsamples = 0;
    bool first = true;
    for (Segment& s : tr.segs) {
        s.build_summary();
        double mn, mx;
        s.minmax(0, s.size(), mn, mx);
        tr.v_min = std::min(tr.v_min, mn);
        tr.v_max = std::max(tr.v_max, mx);
        tr.nsamples += s.size();
        if (first) {
            tr.t_min = s.t0;
            tr.rate = 1.0 / s.dt;
            first = false;
        }
        tr.t_max = s.t_end();
    }
    if (tr.v_min > tr.v_max) tr.v_min = tr.v_max = 0;
    tr.loaded = true;
}

// What one file holds. With data: decoded segments. Without: only the time spans.
struct ScanSeg {
    double t0, t1, rate;
    int64_t n;
};
struct FileChannel {
    Trace meta;
    std::vector<Segment> segs;   // when read with data
    std::vector<ScanSeg> spans;  // when only scanned
};
struct FileContent {
    std::vector<FileChannel> chans;
    std::vector<std::string> messages;
};

// Reads every record of a file and joins records that follow each other straight into
// one segment per channel. Doing the joining here (instead of in libmseed's trace list)
// avoids growing the sample buffer record by record, which is slow for long files.
FileContent read_file_data(const std::string& file) {
    FileContent fc;
    if (!std::ifstream(file).good()) {
        fc.messages.push_back(file + ": file not found or not readable");
        return fc;
    }
    MS3FileParam* fp = nullptr;
    MS3Record* msr = nullptr;
    std::map<std::string, size_t> chan_index;
    std::vector<Segment> cur;  // the segment being extended, for each channel
    std::vector<bool> active;
    int rv, nseg_text = 0;
    while ((rv = ms3_readmsr_r(&fp, &msr, file.c_str(), nullptr, nullptr, MSF_UNPACKDATA, 0)) == MS_NOERROR) {
        if (msr->samprate <= 0 || msr->numsamples <= 0) continue;
        if (msr->sampletype != 'i' && msr->sampletype != 'f' && msr->sampletype != 'd') {
            nseg_text++;
            continue;
        }
        auto it = chan_index.find(msr->sid);
        if (it == chan_index.end()) {
            FileChannel ch;
            parse_sid(msr->sid, ch.meta);
            fc.chans.push_back(std::move(ch));
            cur.emplace_back();
            active.push_back(false);
            it = chan_index.emplace(msr->sid, fc.chans.size() - 1).first;
        }
        const size_t ci = it->second;
        Segment& seg = cur[ci];
        const double t0 = (double)msr->starttime / 1e9, dt = 1.0 / msr->samprate;
        const Segment::Kind kind = msr->sampletype == 'i' ? Segment::Kind::Int32 : Segment::Kind::Float32;
        bool joins = active[ci] && seg.kind == kind && std::fabs(dt - seg.dt) < 1e-6 * dt &&
                     std::fabs(t0 - (seg.t0 + seg.dt * (double)seg.size())) < 0.5 * dt;
        if (!joins) {
            if (active[ci]) fc.chans[ci].segs.push_back(std::move(seg));
            seg = Segment();
            seg.t0 = t0;
            seg.dt = dt;
            seg.kind = kind;
            active[ci] = true;
        }
        const int64_t n = msr->numsamples;
        if (msr->sampletype == 'i') {
            const int32_t* p = (const int32_t*)msr->datasamples;
            seg.i32.insert(seg.i32.end(), p, p + n);
        } else if (msr->sampletype == 'f') {
            const float* p = (const float*)msr->datasamples;
            seg.f32.insert(seg.f32.end(), p, p + n);
        } else {
            const double* p = (const double*)msr->datasamples;
            for (int64_t i = 0; i < n; i++) seg.f32.push_back((float)p[i]);
        }
    }
    if (rv != MS_ENDOFFILE) fc.messages.push_back(file + ": " + ms_errorstr(rv) + " (read what could be read)");
    for (size_t ci = 0; ci < cur.size(); ci++)
        if (active[ci]) fc.chans[ci].segs.push_back(std::move(cur[ci]));
    ms3_readmsr_r(&fp, &msr, nullptr, nullptr, nullptr, 0, 0);  // frees what the reader holds
    if (nseg_text) fc.messages.push_back(file + ": skipped " + std::to_string(nseg_text) + " non-numeric record(s)");
    return fc;
}

FileContent read_file(const std::string& file, bool with_data) {
    if (with_data) return read_file_data(file);
    FileContent fc;
    if (!std::ifstream(file).good()) {
        fc.messages.push_back(file + ": file not found or not readable");
        return fc;
    }
    MS3TraceList* mstl = nullptr;
    int rv = ms3_readtracelist(&mstl, file.c_str(), -1.0, -1.0, 0, with_data ? MSF_UNPACKDATA : 0, 0);
    if (rv != MS_NOERROR && !(mstl && mstl->numtraces > 0)) {
        fc.messages.push_back(file + ": " + ms_errorstr(rv));
        if (mstl) mstl3_free(&mstl, 1);
        return fc;
    }
    if (rv != MS_NOERROR) fc.messages.push_back(file + ": " + ms_errorstr(rv) + " (read what could be read)");
    int nseg_text = 0;
    for (MS3TraceID* id = mstl ? mstl->traces : nullptr; id; id = id->next) {
        FileChannel ch;
        parse_sid(id->sid, ch.meta);
        for (MS3TraceSeg* s = id->first; s; s = s->next) {
            if (s->samprate <= 0) continue;
            if (with_data) {
                if (s->numsamples <= 0) continue;
                Segment seg;
                seg.t0 = (double)s->starttime / 1e9;
                seg.dt = 1.0 / s->samprate;
                if (!convert_samples(s, seg)) {
                    nseg_text++;
                    continue;
                }
                ch.segs.push_back(std::move(seg));
            } else if (s->samplecnt > 0) {
                ch.spans.push_back({(double)s->starttime / 1e9, (double)s->endtime / 1e9, s->samprate, s->samplecnt});
            }
        }
        fc.chans.push_back(std::move(ch));
    }
    if (nseg_text) fc.messages.push_back(file + ": skipped " + std::to_string(nseg_text) + " non-numeric segment(s)");
    if (mstl) mstl3_free(&mstl, 1);
    return fc;
}

}  // namespace

void scan_files(std::vector<Trace>& traces, const std::vector<std::string>& files,
                std::vector<std::string>& messages) {
    std::map<std::string, size_t> index;
    for (size_t i = 0; i < traces.size(); i++) index[traces[i].id] = i;
    int next_order = 0;
    for (const Trace& t : traces) next_order = std::max(next_order, t.file_order + 1);

    std::vector<FileContent> res(files.size());
    parallel_for((int)files.size(), [&](int k) { res[k] = read_file(files[k], false); });

    for (size_t k = 0; k < files.size(); k++) {
        const int order = next_order++;
        for (const std::string& m : res[k].messages) messages.push_back(m);
        for (FileChannel& ch : res[k].chans) {
            if (ch.spans.empty()) continue;
            auto it = index.find(ch.meta.id);
            if (it == index.end()) {
                ch.meta.file = files[k];
                ch.meta.file_order = order;
                ch.meta.t_min = 1e300;
                ch.meta.t_max = -1e300;
                traces.push_back(ch.meta);
                it = index.emplace(ch.meta.id, traces.size() - 1).first;
            }
            Trace& tr = traces[it->second];
            if (std::find(tr.sources.begin(), tr.sources.end(), files[k]) == tr.sources.end()) tr.sources.push_back(files[k]);
            for (const ScanSeg& sp : ch.spans) {
                if (tr.nsamples == 0 && tr.rate == 0) tr.rate = sp.rate;
                tr.t_min = std::min(tr.t_min, sp.t0);
                tr.t_max = std::max(tr.t_max, sp.t1);
                tr.nsamples += sp.n;
            }
        }
    }
    traces.erase(std::remove_if(traces.begin(), traces.end(), [](const Trace& t) { return t.nsamples == 0; }),
                 traces.end());
}

void load_traces(const std::vector<Trace*>& wanted, std::vector<std::string>& messages) {
    std::map<std::string, Trace*> by_id;
    std::set<std::string> fileset;
    for (Trace* t : wanted) {
        if (t->loaded) continue;
        by_id[t->id] = t;
        for (const std::string& f : t->sources) fileset.insert(f);
    }
    if (by_id.empty()) return;
    std::vector<std::string> files(fileset.begin(), fileset.end());
    std::vector<FileContent> res(files.size());
    parallel_for((int)files.size(), [&](int k) { res[k] = read_file(files[k], true); });
    for (FileContent& fc : res) {
        for (const std::string& m : fc.messages) messages.push_back(m);
        for (FileChannel& ch : fc.chans) {
            auto it = by_id.find(ch.meta.id);
            if (it == by_id.end()) continue;  // a channel nobody asked for
            for (Segment& s : ch.segs) it->second->segs.push_back(std::move(s));
        }
    }
    parallel_for((int)by_id.size(), [&](int k) {
        auto it = by_id.begin();
        std::advance(it, k);
        finish_trace(*it->second);
    });
}

void unload_trace(Trace& tr) {
    tr.segs.clear();
    tr.segs.shrink_to_fit();
    tr.loaded = false;
}

LoadResult load_files(const std::vector<std::string>& files) {
    LoadResult r;
    scan_files(r.traces, files, r.messages);
    std::vector<Trace*> all;
    for (Trace& t : r.traces) all.push_back(&t);
    load_traces(all, r.messages);
    return r;
}

std::vector<Segment> extract(const Trace& tr, double t1, double t2) {
    std::vector<Segment> out;
    for (const Segment& s : tr.segs) {
        int64_t i0, i1;
        s.index_range(t1, t2, i0, i1);
        if (i1 <= i0) continue;
        Segment c;
        c.t0 = s.time_at(i0);
        c.dt = s.dt;
        c.kind = s.kind;
        if (s.kind == Segment::Kind::Int32)
            c.i32.assign(s.i32.begin() + i0, s.i32.begin() + i1);
        else
            c.f32.assign(s.f32.begin() + i0, s.f32.begin() + i1);
        out.push_back(std::move(c));
    }
    return out;
}

// ---------------------------------------------------------------------------

std::string format_time(double t, int decimals) {
    double whole = std::floor(t);
    double frac = t - whole;
    int64_t secs = (int64_t)whole;
    // Round the fraction first so 59.9996 does not print as 59.1000.
    double scale = std::pow(10.0, decimals);
    int64_t fr = (int64_t)std::llround(frac * scale);
    if (fr >= (int64_t)scale) {
        fr -= (int64_t)scale;
        secs += 1;
    }
    int64_t days = secs / 86400, rem = secs % 86400;
    if (rem < 0) {
        rem += 86400;
        days -= 1;
    }
    // Days since 1970-01-01 to calendar date (Howard Hinnant's algorithm).
    int64_t z = days + 719468;
    int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    int64_t doe = z - era * 146097;
    int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    int64_t y = yoe + era * 400;
    int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    int64_t mp = (5 * doy + 2) / 153;
    int64_t d = doy - (153 * mp + 2) / 5 + 1;
    int64_t m = mp < 10 ? mp + 3 : mp - 9;
    if (m <= 2) y++;
    char buf[64];
    if (decimals > 0)
        std::snprintf(buf, sizeof buf, "%04lld-%02lld-%02lld %02lld:%02lld:%02lld.%0*lld", (long long)y,
                      (long long)m, (long long)d, (long long)(rem / 3600), (long long)(rem / 60 % 60),
                      (long long)(rem % 60), decimals, (long long)fr);
    else
        std::snprintf(buf, sizeof buf, "%04lld-%02lld-%02lld %02lld:%02lld:%02lld", (long long)y, (long long)m,
                      (long long)d, (long long)(rem / 3600), (long long)(rem / 60 % 60), (long long)(rem % 60));
    return buf;
}

}  // namespace msv
