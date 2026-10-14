#include "app.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "imgui.h"
#include "implot.h"
#include "parallel.h"
#include "stats.h"

namespace msv {

namespace {

// Extra height the bottom plot needs for its time/frequency labels.
}  // namespace
float axis_label_height() { return ImGui::GetFontSize() * 2.6f + 12.0f; }

namespace {
const char* kDetrendNames[] = {"none", "remove mean", "remove linear trend"};

}  // namespace

// ---------------------------------------------------------------------------
// Setup and bookkeeping

App::App(const std::vector<std::string>& files) {
    load_config(cfg);
    load_filters(filters_);
    load_calib(calib_);
    proc_.detrend = cfg.detrend;
    proc_.reverse = cfg.reverse_polarity;
    mag_.overlay = cfg.overlay_magnify;
    split_mag_.overlay = cfg.split_overlay_magnify;
    for (SpecState* s : {&spec_, &split_spec_}) {
        s->welch = cfg.welch;
        s->nwin = cfg.welch_window;
        s->log_x = cfg.spec_log_x;
        s->log_y = cfg.spec_log_y;
        s->normalize = cfg.spec_normalize;
    }
    spec_.overlay = cfg.overlay_spectra;
    split_spec_.overlay = cfg.split_overlay_spectra;

    ImPlotInputMap& m = ImPlot::GetInputMap();
    m.Pan = ImGuiMouseButton_Right;    // left button is for picking times
    m.Fit = ImGuiMouseButton_Middle;
    m.Select = ImGuiMouseButton_Middle;
    m.Menu = ImGuiMouseButton_Middle;
    ImPlotStyle& ps = ImPlot::GetStyle();
    ps.UseISO8601 = true;
    ps.Use24HourClock = true;
    ps.PlotPadding = ImVec2(4, 1);  // little space around plots, so stacked plots use the screen well

    add_files(files);
}

void App::save_all() {
    cfg.detrend = proc_.detrend;
    cfg.reverse_polarity = proc_.reverse;
    cfg.overlay_magnify = mag_.overlay;
    cfg.overlay_spectra = spec_.overlay;
    cfg.split_overlay_magnify = split_mag_.overlay;
    cfg.split_overlay_spectra = split_spec_.overlay;
    cfg.welch = spec_.welch;
    cfg.welch_window = spec_.nwin;
    cfg.spec_log_x = spec_.log_x;
    cfg.spec_log_y = spec_.log_y;
    cfg.spec_normalize = spec_.normalize;
    save_config(cfg);
    save_filters(filters_);
    save_calib(calib_);
}

void App::add_files(const std::vector<std::string>& files) {
    if (files.empty()) return;
    n_files_ += (int)files.size();
    ::msv::scan_files(traces_, files, messages_);
    topts_.resize(traces_.size());
    calib_sel_.assign(traces_.size(), 0);
    last_used_.assign(traces_.size(), 0);
    rebuild_order();
    reset_views();
    clear_caches();
}

void App::rebuild_order() {
    order_ = sort_order(traces_, cfg.sort);
    page_ = std::max(0, std::min(page_, n_pages() - 1));
}

void App::reset_views() {
    double lo, hi;
    data_range(lo, hi);
    main_lim_[0] = lo;
    main_lim_[1] = hi;
    sel_state_ = 0;
    have_win_ = false;
    pick_state_ = 0;
    sel_applied_ = false;
}

void App::clear_caches() {
    proc_cache_.clear();
    spec_cache_.clear();
    spec_stale_.clear();
    rms_cache_.clear();
    stats_cache_.clear();
    stats_stale_.clear();
}

double App::off(int ti) const { return cfg.relative_time ? traces_[ti].t_min : 0.0; }

void App::data_range(double& lo, double& hi) const {
    if (traces_.empty()) {
        lo = 0;
        hi = 1;
        return;
    }
    if (cfg.relative_time) {
        lo = 0;
        hi = 0;
        for (const Trace& t : traces_) hi = std::max(hi, t.t_max - t.t_min);
    } else {
        lo = traces_[0].t_min;
        hi = traces_[0].t_max;
        for (const Trace& t : traces_) {
            lo = std::min(lo, t.t_min);
            hi = std::max(hi, t.t_max);
        }
    }
    if (hi <= lo) hi = lo + 1;
}

std::string App::fmt_disp(double d) const {
    if (cfg.relative_time) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.3f s", d);
        return buf;
    }
    return format_time(d) + " UTC";
}

int App::n_pages() const {
    int n = std::max(1, cfg.n_per_screen);
    return std::max(1, ((int)order_.size() + n - 1) / n);
}

std::vector<int> App::page_traces(bool only_on) const {
    std::vector<int> r;
    int n = std::max(1, cfg.n_per_screen);
    for (int k = page_ * n; k < std::min((int)order_.size(), (page_ + 1) * n); k++) {
        int ti = order_[k];
        if (!only_on || traces_[ti].on) r.push_back(ti);
    }
    return r;
}

const ProcOpts& App::opts_for(int ti) const { return topts_[ti].custom ? topts_[ti].o : proc_; }

ImVec4 App::trace_color(int ti) const {
    uint32_t c = cfg.colors[(size_t)ti % kNumTraceColors];
    return ImVec4(((c >> 16) & 255) / 255.0f, ((c >> 8) & 255) / 255.0f, (c & 255) / 255.0f, 1.0f);
}

ImVec4 App::mark_color() const {
    uint32_t c = cfg.colors[kNumTraceColors];
    return ImVec4(((c >> 16) & 255) / 255.0f, ((c >> 8) & 255) / 255.0f, (c & 255) / 255.0f, 1.0f);
}

// A list of powers of two to pick from (so that a stray click cannot make it 1025).
bool App::combo_pow2(const char* label, int& v, int lo, int hi) {
    bool changed = false;
    v = std::max(lo, std::min(hi, snap_pow2(v)));
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
    if (ImGui::BeginCombo(label, std::to_string(v).c_str())) {
        for (int p = lo; p <= hi; p *= 2)
            if (ImGui::Selectable(std::to_string(p).c_str(), p == v)) {
                v = p;
                changed = true;
            }
        ImGui::EndCombo();
    }
    return changed;
}

const Calib& App::calib_for(int ti) const {
    static const Calib none;
    auto it = calib_.find(traces_[ti].id);
    return it == calib_.end() ? none : it->second;
}

std::string App::units_of(const std::vector<int>& rows) const {
    std::string u;
    for (int ti : rows) {
        const std::string& x = calib_for(ti).units;
        if (u.empty()) u = x;
        else if (u != x) return "mixed units";
    }
    return u.empty() ? "counts" : u;
}

bool App::busy() const {
    // Spectra wait 0.15 s after the view stops moving, then appear; keep redrawing until then.
    return loading_pending_ || (have_win_ && ImGui::GetTime() - view_changed_ < 0.4);
}

bool App::view_stable() const { return ImGui::GetTime() - view_changed_ > 0.15; }

// ---------------------------------------------------------------------------
// Processing

std::string App::proc_key(int ti, bool wf) const {
    const ProcOpts& o = opts_for(ti);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%d|%.6f|%.6f|%d|%d|%d|%.9g|", ti, off(ti) + win_[0], off(ti) + win_[1], o.detrend,
                  (int)o.reverse, (int)wf, calib_for(ti).scale);
    std::string k = buf;
    if (wf && o.filter >= 0 && o.filter < (int)filters_.size()) {
        const FilterSpec& f = filters_[o.filter];
        std::snprintf(buf, sizeof buf, "%d,%d,%.9g,%.9g,%d", (int)f.type, f.order, f.f1, f.f2, (int)f.zero_phase);
        k += buf;
    }
    return k;
}

std::shared_ptr<const Processed> App::compute_processed(int ti, bool wf) const {
    auto P = std::make_shared<Processed>();
    const Trace& tr = traces_[ti];
    const ProcOpts& o = opts_for(ti);
    const FilterSpec* f = nullptr;
    if (wf && o.filter >= 0 && o.filter < (int)filters_.size()) {
        f = &filters_[o.filter];
        std::string why = check_filter(*f, tr.rate);
        if (!why.empty()) {
            P->warning = "filter not applied: " + why;
            f = nullptr;
        }
    }
    for (const Segment& raw : extract(tr, off(ti) + win_[0], off(ti) + win_[1])) {
        std::vector<double> v = raw.doubles();
        const double scale = calib_for(ti).scale;
        if (scale != 1.0)
            for (double& d : v) d *= scale;
        detrend(v, (Detrend)o.detrend);
        if (f) apply_filter(*f, 1.0 / raw.dt, v);
        if (o.reverse)
            for (double& d : v) d = -d;
        Segment s;
        s.t0 = raw.t0;
        s.dt = raw.dt;
        s.set(v);
        s.build_summary();
        P->segs.push_back(std::move(s));
    }
    return P;
}

void App::ensure_processed(const std::vector<int>& rows, bool wf) {
    if (!have_win_) return;
    std::vector<int> missing;
    for (int ti : rows)
        if (!proc_cache_.count(proc_key(ti, wf))) missing.push_back(ti);
    if (missing.empty()) return;
    if (proc_cache_.size() + missing.size() > (size_t)(3 * std::max(1, cfg.n_per_screen) + 12)) {
        proc_cache_.clear();
        spec_cache_.clear();
        spec_stale_.clear();
        rms_cache_.clear();
        missing = rows;
    }
    std::vector<std::shared_ptr<const Processed>> res(missing.size());
    parallel_for((int)missing.size(), [&](int k) { res[k] = compute_processed(missing[k], wf); });
    for (size_t k = 0; k < missing.size(); k++) proc_cache_[proc_key(missing[k], wf)] = res[k];
}

std::shared_ptr<const Processed> App::proc_for(int ti, bool wf) const {
    auto it = proc_cache_.find(proc_key(ti, wf));
    return it == proc_cache_.end() ? nullptr : it->second;
}

double App::rms_for(int ti, bool wf) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "|%.6f|%.6f", view_[0], view_[1]);
    std::string key = proc_key(ti, wf) + buf;
    auto it = rms_cache_.find(key);
    if (it != rms_cache_.end()) return it->second;
    auto P = proc_for(ti, wf);
    if (!P) return std::nan("");
    double sum = 0;
    int64_t n = 0;
    for (const Segment& s : P->segs) {
        int64_t i0, i1;
        s.index_range(off(ti) + view_[0], off(ti) + view_[1], i0, i1);
        for (int64_t i = i0; i < i1; i++) {
            double v = s.at(i);
            sum += v * v;
        }
        n += i1 - i0;
    }
    double r = n ? std::sqrt(sum / (double)n) : std::nan("");
    if (rms_cache_.size() > 500) rms_cache_.clear();
    rms_cache_[key] = r;
    return r;
}

namespace {

// Thins a spectrum for drawing: neighbouring points closer than 0.5% in
// frequency are averaged together.
SpecResult make_spec_result(const Spectrum& sp) {
    SpecResult r;
    const double kStep = 1.005;
    size_t i = 0;
    while (i < sp.f.size()) {
        double f0 = sp.f[i], fsum = 0, psum = 0;
        size_t cnt = 0;
        while (i < sp.f.size() && (cnt == 0 || sp.f[i] < f0 * kStep)) {
            fsum += sp.f[i];
            psum += sp.p[i];
            cnt++;
            i++;
        }
        r.pts.x.push_back(fsum / (double)cnt);
        r.pts.y.push_back(std::max(psum / (double)cnt, 1e-30));
        r.pmax = std::max(r.pmax, r.pts.y.back());
    }
    if (r.pmax < 1e-20) {  // no signal at all (e.g. a constant trace): nothing worth drawing
        r.pts.clear();
        r.pmax = 0;
    }
    return r;
}

}  // namespace

std::vector<std::shared_ptr<const StatsResult>> App::stats_for(const std::vector<int>& rows, const MagState& ms,
                                                               int nbins, const char* panel) {
    std::vector<std::shared_ptr<const StatsResult>> out(rows.size());
    const bool stable = view_stable();
    std::vector<std::string> keys(rows.size()), stale_keys(rows.size());
    std::vector<int> todo;
    for (size_t r = 0; r < rows.size(); r++) {
        int ti = rows[r];
        char buf[160];
        std::snprintf(buf, sizeof buf, "|v%.6f|%.6f|w%.6g|%d", off(ti) + view_[0], off(ti) + view_[1], ms.st_win, nbins);
        keys[r] = proc_key(ti, true) + buf;
        stale_keys[r] = std::string(panel) + "|" + std::to_string(ti);
        auto it = stats_cache_.find(keys[r]);
        if (it != stats_cache_.end()) {
            out[r] = it->second;
            stats_stale_[stale_keys[r]] = it->second;
        } else if (stable) {
            todo.push_back((int)r);
        } else {
            auto st = stats_stale_.find(stale_keys[r]);
            if (st != stats_stale_.end()) out[r] = st->second;
        }
    }
    if (todo.empty()) return out;
    if (stats_cache_.size() > 300) stats_cache_.clear();
    std::vector<std::shared_ptr<const StatsResult>> res(todo.size());
    parallel_for((int)todo.size(), [&](int k) {
        int ti = rows[todo[k]];
        auto P = proc_for(ti, true);
        if (!P) return;
        auto R = std::make_shared<StatsResult>();
        for (const Segment& seg : P->segs) {
            MovingStats m = moving_stats(seg, off(ti) + view_[0], off(ti) + view_[1], ms.st_win, nbins);
            Series a, b, c;
            for (size_t i = 0; i < m.t.size(); i++) {
                a.x.push_back(m.t[i]);
                a.y.push_back(m.rms[i]);
                if (!std::isnan(m.skew[i])) {
                    b.x.push_back(m.t[i]);
                    b.y.push_back(m.skew[i]);
                    c.x.push_back(m.t[i]);
                    c.y.push_back(m.kurt[i]);
                }
            }
            R->rms.push_back(std::move(a));
            R->skew.push_back(std::move(b));
            R->kurt.push_back(std::move(c));
        }
        res[k] = R;
    });
    for (size_t k = 0; k < todo.size(); k++) {
        if (!res[k]) continue;
        stats_cache_[keys[todo[k]]] = res[k];
        stats_stale_[stale_keys[todo[k]]] = res[k];
        out[todo[k]] = res[k];
    }
    return out;
}

std::vector<std::shared_ptr<const SpecResult>> App::spectra_for(const std::vector<int>& rows, bool wf,
                                                                const SpecState& ss, const char* panel) {
    std::vector<std::shared_ptr<const SpecResult>> out(rows.size());
    const bool stable = view_stable();
    std::vector<std::string> keys(rows.size()), stale_keys(rows.size());
    std::vector<int> todo;
    for (size_t r = 0; r < rows.size(); r++) {
        int ti = rows[r];
        char buf[160];
        std::snprintf(buf, sizeof buf, "|v%.6f|%.6f|%d|%d", off(ti) + view_[0], off(ti) + view_[1], (int)ss.welch,
                      ss.nwin);
        keys[r] = proc_key(ti, wf) + buf;
        stale_keys[r] = std::string(panel) + "|" + std::to_string(ti) + "|" + std::to_string((int)wf);
        auto it = spec_cache_.find(keys[r]);
        if (it != spec_cache_.end()) {
            out[r] = it->second;
            spec_stale_[stale_keys[r]] = it->second;
        } else if (stable) {
            todo.push_back((int)r);
        } else {
            auto st = spec_stale_.find(stale_keys[r]);
            if (st != spec_stale_.end()) out[r] = st->second;
        }
    }
    if (todo.empty()) return out;
    if (spec_cache_.size() > 300) spec_cache_.clear();
    std::vector<std::shared_ptr<const SpecResult>> res(todo.size());
    parallel_for((int)todo.size(), [&](int k) {
        int ti = rows[todo[k]];
        auto P = proc_for(ti, wf);
        if (!P) return;
        // Use the longest gap-free piece inside the view.
        const Segment* best = nullptr;
        int64_t b0 = 0, b1 = 0;
        for (const Segment& s : P->segs) {
            int64_t i0, i1;
            s.index_range(off(ti) + view_[0], off(ti) + view_[1], i0, i1);
            if (i1 - i0 > b1 - b0) {
                best = &s;
                b0 = i0;
                b1 = i1;
            }
        }
        auto R = std::make_shared<SpecResult>();
        if (best && b1 - b0 >= 8) {
            std::vector<double> v = best->doubles(b0, b1);
            Spectrum sp = power_spectrum(v.data(), v.size(), 1.0 / best->dt, ss.welch, ss.nwin);
            *R = make_spec_result(sp);
        }
        res[k] = R;
    });
    for (size_t k = 0; k < todo.size(); k++) {
        if (!res[k]) continue;
        spec_cache_[keys[todo[k]]] = res[k];
        spec_stale_[stale_keys[todo[k]]] = res[k];
        out[todo[k]] = res[k];
    }
    return out;
}

// ---------------------------------------------------------------------------
// Frame

// Reads the samples of the traces on the current page if they are not in memory yet.
// Returns false (after drawing a message) on the frame before the reading starts, so the
// message is visible while the program is busy. Old pages are dropped to limit memory.
bool App::ensure_page_loaded() {
    last_used_.resize(traces_.size(), 0);
    std::vector<Trace*> need;
    int64_t samples = 0;
    for (int ti : page_traces(false)) {
        last_used_[ti] = frame_no_;
        if (!traces_[ti].loaded) {
            need.push_back(&traces_[ti]);
            samples += traces_[ti].nsamples;
        }
    }
    if (need.empty()) {
        loading_pending_ = false;
        return true;
    }
    if (!loading_pending_) {
        loading_pending_ = true;
        ImGui::Spacing();
        ImGui::Text("Reading %d channel(s), %.1f million samples ...", (int)need.size(), samples / 1e6);
        return false;
    }
    loading_pending_ = false;
    ::msv::load_traces(need, messages_);
    // Keep at most two pages of traces in memory: drop the ones used longest ago.
    const size_t keep = 2 * (size_t)std::max(1, cfg.n_per_screen);
    std::vector<int> loaded;
    for (size_t i = 0; i < traces_.size(); i++)
        if (traces_[i].loaded) loaded.push_back((int)i);
    if (loaded.size() > keep) {
        std::sort(loaded.begin(), loaded.end(), [&](int a, int b) { return last_used_[a] < last_used_[b]; });
        for (size_t k = 0; k + keep < loaded.size(); k++)
            if (last_used_[loaded[k]] != frame_no_) unload_trace(traces_[loaded[k]]);
    }
    return true;
}

void App::frame() {
    frame_no_++;
    // Note when the magnify axes move, so spectra can wait until they settle.
    if (view_[0] != last_view_[0] || view_[1] != last_view_[1]) {
        last_view_[0] = view_[0];
        last_view_[1] = view_[1];
        view_changed_ = ImGui::GetTime();
    }
    grow_window_if_needed();
    draw_top_bar();
    if (show_config_) {
        draw_config();
        goto_tab_ = -1;
        return;
    }
    if (traces_.empty()) {
        ImGui::Spacing();
        ImGui::TextWrapped(
            "No data loaded. Drop miniSEED files onto this window, or start the program with file names:\n"
            "    miniseed_viewer file1.mseed file2.mseed ...");
        return;
    }
    if (!ensure_page_loaded()) return;
    if (ImGui::BeginTabBar("tabs")) {
        auto flag = [&](int i) { return goto_tab_ == i ? ImGuiTabItemFlags_SetSelected : 0; };
        if (ImGui::BeginTabItem("Main", nullptr, flag(0))) {
            draw_main();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Magnify", nullptr, flag(1))) {
            apply_main_selection();
            draw_magnify("##mag", ImGui::GetContentRegionAvail().y, mag_);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Spectra", nullptr, flag(2))) {
            apply_main_selection();
            draw_spectra("##spec", ImGui::GetContentRegionAvail().y, spec_);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Spectrogram", nullptr, flag(3))) {
            apply_main_selection();
            draw_specgram("##sg", ImGui::GetContentRegionAvail().y, sg_);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Split", nullptr, flag(4))) {
            apply_main_selection();
            draw_split();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    goto_tab_ = -1;
}

void App::draw_top_bar() {
    ImGui::Text("%d file(s), %d channel(s)", n_files_, (int)traces_.size());
    if (n_pages() > 1 && !show_config_) {
        ImGui::SameLine(0, 30);
        ImGui::BeginDisabled(page_ == 0);
        if (ImGui::Button("< Previous")) page_ = std::max(0, page_ - 1);
        ImGui::EndDisabled();
        ImGui::SameLine();
        const int n = std::max(1, cfg.n_per_screen);
        ImGui::Text("page %d of %d (channels %d-%d of %d)", page_ + 1, n_pages(), page_ * n + 1,
                    std::min((int)order_.size(), (page_ + 1) * n), (int)order_.size());
        ImGui::SameLine();
        ImGui::BeginDisabled(page_ + 1 >= n_pages());
        if (ImGui::Button("Next >")) page_ = std::min(n_pages() - 1, page_ + 1);
        ImGui::EndDisabled();
        if (!ImGui::GetIO().WantTextInput) {
            if (ImGui::IsKeyPressed(ImGuiKey_PageUp)) page_ = std::max(0, page_ - 1);
            if (ImGui::IsKeyPressed(ImGuiKey_PageDown)) page_ = std::min(n_pages() - 1, page_ + 1);
        }
    }
    const char* label = show_config_ ? "Back" : "Settings";
    ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - ImGui::CalcTextSize(label).x -
                    2 * ImGui::GetStyle().FramePadding.x);
    if (ImGui::Button(label)) {
        show_config_ = !show_config_;
        if (!show_config_) save_all();
    }
    if (!messages_.empty()) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.75f, 0.3f, 1));
        for (const std::string& m : messages_) ImGui::TextUnformatted(m.c_str());
        ImGui::PopStyleColor();
        if (ImGui::SmallButton("Dismiss messages")) messages_.clear();
    }
    ImGui::Separator();
}

void App::need_window_message() {
    ImGui::Spacing();
    ImGui::TextWrapped(
        "No time window chosen yet. On the Main tab, click once on a plot to set t1 and again to set t2; then this tab "
        "will show that window.");
    if (ImGui::Button("Use the full time range")) {
        double lo, hi;
        data_range(lo, hi);
        apply_window(lo, hi);
    }
}

// Sets the magnify view. The processed window gets some margin around it, so that
// filter start-up effects stay out of sight and small pans need no recalculation.
void App::apply_window(double a, double b) {
    view_[0] = a;
    view_[1] = b;
    double lo, hi;
    data_range(lo, hi);
    double span = b - a;
    win_[0] = std::min(a, std::max(lo, a - 0.5 * span));
    win_[1] = std::max(b, std::min(hi, b + 0.5 * span));
    have_win_ = true;
    pick_state_ = 0;
}

// Called when a tab that needs a time window is shown: if the selection made on the Main tab is new,
// the magnify view is set to it (so Magnify, Spectra, Spectrogram and Split all follow Main).
void App::apply_main_selection() {
    if (sel_state_ != 2) return;
    if (sel_applied_ && sel_[0] == applied_sel_[0] && sel_[1] == applied_sel_[1]) return;
    apply_window(sel_[0], sel_[1]);
    applied_sel_[0] = sel_[0];
    applied_sel_[1] = sel_[1];
    sel_applied_ = true;
}

// When the view has moved outside the processed window (zoomed or panned out),
// process the new region once the view has stopped moving.
void App::grow_window_if_needed() {
    if (!have_win_ || !view_stable()) return;
    if (view_[0] >= win_[0] - 1e-9 && view_[1] <= win_[1] + 1e-9) return;
    apply_window(view_[0], view_[1]);
}

// ---------------------------------------------------------------------------
// Shared small widgets

void App::combo_filter(const char* label, int& idx) {
    const char* cur = (idx >= 0 && idx < (int)filters_.size()) ? filters_[idx].name.c_str() : "none";
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
    if (ImGui::BeginCombo(label, cur)) {
        if (ImGui::Selectable("none", idx < 0)) idx = -1;
        for (int i = 0; i < (int)filters_.size(); i++)
            if (ImGui::Selectable((filters_[i].name + "##f" + std::to_string(i)).c_str(), idx == i)) idx = i;
        ImGui::EndCombo();
    }
}

void App::combo_detrend(const char* label, int& d) {
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
    ImGui::Combo(label, &d, kDetrendNames, 3);
}

void App::trace_label(int ti, bool checkbox, bool settings_button, const char* line2) {
    Trace& tr = traces_[ti];
    ImVec4 col = trace_color(ti);
    if (!tr.on) col = ImVec4(0.5f, 0.5f, 0.5f, 1);
    if (checkbox) {
        // The name is the check box's label, so clicking the name works like clicking the box.
        ImGui::PushStyleColor(ImGuiCol_Text, col);
        ImGui::Checkbox(tr.id.c_str(), &tr.on);
        ImGui::PopStyleColor();
    } else {
        ImGui::TextColored(col, "%s", tr.id.c_str());
    }
    if (line2) ImGui::TextUnformatted(line2);
    if (settings_button) {
        TraceOpts& to = topts_[ti];
        std::string lbl = std::string(to.custom ? "custom*" : "settings") + "##s";
        if (ImGui::SmallButton(lbl.c_str())) ImGui::OpenPopup("proc");
        if (ImGui::BeginPopup("proc")) {
            ImGui::TextUnformatted(tr.id.c_str());
            ImGui::Checkbox("Use different settings for this trace", &to.custom);
            if (to.custom) {
                combo_filter("filter", to.o.filter);
                combo_detrend("detrend", to.o.detrend);
                ImGui::Checkbox("reverse polarity", &to.o.reverse);
                if (ImGui::Button("Copy the global settings here")) to.o = proc_;
            }
            ImGui::EndPopup();
        }
    }
}

// ---------------------------------------------------------------------------
// Main tab

void App::main_click(double x) {
    if (sel_state_ != 1) {
        sel_[0] = x;
        sel_state_ = 1;
    } else {
        sel_[1] = x;
        if (sel_[1] < sel_[0]) std::swap(sel_[0], sel_[1]);
        sel_state_ = 2;
    }
}

void App::draw_main() {
    if (ImGui::Button("Reset view")) {
        double lo, hi;
        data_range(lo, hi);
        main_lim_[0] = lo;
        main_lim_[1] = hi;
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear selection")) sel_state_ = 0;
    ImGui::SameLine();
    if (ImGui::Button("Calibration...")) ImGui::OpenPopup("Calibration");
    calib_editor();
    ImGui::SameLine(0, 30);
    if (sel_state_ == 0)
        ImGui::TextUnformatted("Click on a plot to set t1, click again to set t2. Then open any other tab.");
    else if (sel_state_ == 1)
        ImGui::Text("t1 = %s   (click again to set t2)", fmt_disp(sel_[0]).c_str());
    else
        ImGui::Text("t1 = %s    t2 = %s    (%.3f s)   Open the Magnify, Spectra, Spectrogram or Split tab to see it.",
                    fmt_disp(sel_[0]).c_str(), fmt_disp(sel_[1]).c_str(), sel_[1] - sel_[0]);

    ImGui::TextUnformatted("Scale:");
    ImGui::SameLine();
    if (ImGui::RadioButton("each trace", !cfg.main_global_norm)) cfg.main_global_norm = false;
    ImGui::SameLine();
    if (ImGui::RadioButton("all traces together", cfg.main_global_norm)) cfg.main_global_norm = true;
    ImGui::SameLine();
    ImGui::Checkbox("fit to the part on screen", &cfg.main_norm_visible);
    ImGui::SameLine(0, 30);
    ImGui::TextUnformatted("Time axis:");
    ImGui::SameLine();
    bool changed = false;
    if (ImGui::RadioButton("clock time (UTC)", !cfg.relative_time) && cfg.relative_time) {
        cfg.relative_time = false;
        changed = true;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("seconds since start of trace", cfg.relative_time) && !cfg.relative_time) {
        cfg.relative_time = true;
        changed = true;
    }
    if (changed) {
        reset_views();
        clear_caches();
    }

    // ---- gather what to draw
    std::vector<int> rows = page_traces(false);
    const float label_w = ImGui::GetFontSize() * 11.0f;
    const float row_h =
        std::max(30.0f, std::floor((ImGui::GetContentRegionAvail().y - axis_label_height()) / std::max<size_t>(1, rows.size())));
    const int nbins = std::max(100, std::min(4000, (int)(ImGui::GetContentRegionAvail().x - label_w - 20)));

    struct Row {
        std::vector<Series> s;
        double lo = 0, hi = 1;
    };
    std::vector<Row> data(rows.size());
    std::map<std::string, std::pair<double, double>> common;  // shared range for each kind of units
    for (size_t r = 0; r < rows.size(); r++) {
        const Trace& tr = traces_[rows[r]];
        const double o = off(rows[r]);
        const double scale = calib_for(rows[r]).scale;  // physical units per count
        Row& d = data[r];
        double vlo = 1e300, vhi = -1e300;
        for (const Segment& seg : tr.segs) {
            Series ser;
            seg.envelope(o + main_lim_[0], o + main_lim_[1], nbins, ser);
            if (ser.x.empty()) continue;
            for (size_t i = 0; i < ser.x.size(); i++) {
                ser.y[i] *= scale;
                vlo = std::min(vlo, ser.y[i]);
                vhi = std::max(vhi, ser.y[i]);
                ser.x[i] -= o;
            }
            d.s.push_back(std::move(ser));
        }
        if (cfg.main_norm_visible && vlo <= vhi) {
            d.lo = vlo;
            d.hi = vhi;
        } else {
            d.lo = std::min(tr.v_min * scale, tr.v_max * scale);
            d.hi = std::max(tr.v_min * scale, tr.v_max * scale);
        }
        if (tr.on) {  // traces that are switched off do not take part in the common scale
            auto it = common.emplace(calib_for(rows[r]).units, std::make_pair(1e300, -1e300)).first;
            it->second.first = std::min(it->second.first, d.lo);
            it->second.second = std::max(it->second.second, d.hi);
        }
    }
    for (size_t r = 0; r < rows.size(); r++) {
        Row& d = data[r];
        if (cfg.main_global_norm && traces_[rows[r]].on) {
            const auto& c = common[calib_for(rows[r]).units];
            d.lo = c.first;
            d.hi = c.second;
        }
        double span = d.hi - d.lo;
        double big = std::max(std::fabs(d.lo), std::fabs(d.hi));
        if (span <= 1e-12 * big || big == 0) {  // flat trace: give it some room
            double m = big > 0 ? big : 1.0;
            d.lo -= m;
            d.hi += m;
        } else {
            d.lo -= 0.05 * span;
            d.hi += 0.05 * span;
        }
    }

    // ---- draw
    bool clicked = false;
    double click_x = 0;
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
    bool aligned = ImPlot::BeginAlignedPlots("##align_main");
    for (size_t r = 0; r < rows.size(); r++) {
        int ti = rows[r];
        const Trace& tr = traces_[ti];
        const std::string& units = calib_for(ti).units;
        ImGui::PushID(ti);
        ImGui::BeginChild("lbl", ImVec2(label_w, row_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
        trace_label(ti, true, false, units == "counts" ? nullptr : units.c_str());
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        const bool last = (r + 1 == rows.size());
        if (ImPlot::BeginPlot("##p", ImVec2(-1, last ? row_h + axis_label_height() : row_h),
                              ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect |
                                  ImPlotFlags_NoFrame)) {
            ImPlot::SetupAxis(ImAxis_X1, nullptr,
                              ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight |
                                  (last ? 0 : ImPlotAxisFlags_NoTickLabels));
            ImPlot::SetupAxis(ImAxis_Y1, nullptr,
                              ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoGridLines |
                                  ImPlotAxisFlags_NoHighlight | ImPlotAxisFlags_Lock);
            ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
            ImPlot::SetupAxisLinks(ImAxis_X1, &main_lim_[0], &main_lim_[1]);
            ImPlot::SetupAxisLimits(ImAxis_Y1, data[r].lo, data[r].hi, ImPlotCond_Always);
            ImPlotSpec spec;
            spec.LineColor = trace_color(ti);
            if (!tr.on) spec.LineColor.w = 0.3f;  // switched-off traces stay visible, but faded
            spec.Flags = ImPlotItemFlags_NoLegend;
            for (const Series& s : data[r].s) ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);

            ImPlotSpec mark;
            mark.LineColor = mark_color();
            mark.LineWeight = 1.5f;
            mark.Flags = ImPlotItemFlags_NoLegend;
            if (sel_state_ >= 1) ImPlot::PlotInfLines("##t1", &sel_[0], 1, mark);
            if (sel_state_ == 2) {
                ImPlot::PlotInfLines("##t2", &sel_[1], 1, mark);
                ImPlotRect lim = ImPlot::GetPlotLimits();
                ImVec2 p0 = ImPlot::PlotToPixels(sel_[0], lim.Y.Max), p1 = ImPlot::PlotToPixels(sel_[1], lim.Y.Min);
                ImPlot::PushPlotClipRect();
                ImVec4 mc = mark_color();
                mc.w = 0.16f;
                ImPlot::GetPlotDrawList()->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32(mc));
                ImPlot::PopPlotClipRect();
            }
            if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                ImVec2 d = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
                if (d.x == 0 && d.y == 0) {
                    clicked = true;
                    click_x = ImPlot::GetPlotMousePos().x;
                }
            }
            ImPlot::EndPlot();
        }
        ImGui::PopID();
    }
    if (aligned) ImPlot::EndAlignedPlots();
    ImGui::PopStyleVar();
    if (clicked) main_click(click_x);
}

// ---------------------------------------------------------------------------
// Settings screen

void App::draw_config() {
    ImGui::BeginChild("##config", ImVec2(0, 0));
    ImGui::SeparatorText("Layout");
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
    if (ImGui::InputInt("Traces per screen", &cfg.n_per_screen)) {
        cfg.n_per_screen = std::max(1, std::min(64, cfg.n_per_screen));
        page_ = std::max(0, std::min(page_, n_pages() - 1));
    }
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 10);
    ImGui::SliderFloat("Text size", &cfg.ui_scale, 0.75f, 2.5f, "%.2fx");

    ImGui::SeparatorText("Main tab");
    ImGui::Checkbox("Scale all traces together (instead of each trace separately)", &cfg.main_global_norm);
    ImGui::Checkbox("Fit the scale to the part on screen (instead of the whole trace)", &cfg.main_norm_visible);
    bool rel = cfg.relative_time;
    if (ImGui::Checkbox("Time axis: seconds since start of trace (instead of clock time)", &rel)) {
        cfg.relative_time = rel;
        reset_views();
        clear_caches();
    }

    ImGui::SeparatorText("Starting values for magnify and spectra tabs");
    combo_detrend("Detrend", proc_.detrend);
    ImGui::Checkbox("Reverse polarity", &proc_.reverse);
    ImGui::Checkbox("Spectra: Welch average (instead of one spectrum)", &spec_.welch);
    combo_pow2("Spectra: Welch window (samples)", spec_.nwin, 64, 65536);
    ImGui::Checkbox("Spectra: log frequency axis", &spec_.log_x);
    ImGui::Checkbox("Spectra: log power axis", &spec_.log_y);
    ImGui::Checkbox("Spectra: normalize", &spec_.normalize);
    ImGui::TextDisabled("(Changes here apply to the Spectra tab now and to the Split tab next time the program starts.)");

    ImGui::SeparatorText("Colors");
    ImGui::TextWrapped(
        "Trace colors are used in turn for the channels. The defaults are chosen to be easy to tell apart with the "
        "common kinds of color blindness. Click a swatch to change it.");
    {
        std::vector<uint32_t> defaults = default_colors();
        for (int i = 0; i <= kNumTraceColors; i++) {
            uint32_t& c = cfg.colors[i];
            float rgb[3] = {((c >> 16) & 255) / 255.0f, ((c >> 8) & 255) / 255.0f, (c & 255) / 255.0f};
            ImGui::PushID(i);
            if (ImGui::ColorEdit3(i < kNumTraceColors ? ("Trace " + std::to_string(i + 1)).c_str() : "t1/t2 marks", rgb,
                                  ImGuiColorEditFlags_NoInputs))
                c = ((uint32_t)std::lround(rgb[0] * 255) << 16) | ((uint32_t)std::lround(rgb[1] * 255) << 8) |
                    (uint32_t)std::lround(rgb[2] * 255);
            ImGui::PopID();
            if (i < kNumTraceColors && i != 3 && i != 7) ImGui::SameLine(0, 20);
        }
        if (ImGui::Button("Restore original colors")) cfg.colors = defaults;
    }

    ImGui::SeparatorText("Order of channels");
    ImGui::TextWrapped(
        "Channels are sorted by the checked items below. The top item is the most important; later items break ties.");
    bool sort_changed = false;
    for (int i = 0; i < (int)cfg.sort.size(); i++) {
        SortItem& s = cfg.sort[i];
        ImGui::PushID(i);
        sort_changed |= ImGui::Checkbox("##use", &s.use);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
        ImGui::TextUnformatted(sort_key_name(s.key));
        ImGui::SameLine(ImGui::GetFontSize() * 17);
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
        int dir = s.descending ? 1 : 0;
        if (ImGui::Combo("##dir", &dir, "ascending\0descending\0")) {
            s.descending = dir == 1;
            sort_changed = true;
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("##up", ImGuiDir_Up) && i > 0) {
            std::swap(cfg.sort[i], cfg.sort[i - 1]);
            sort_changed = true;
        }
        ImGui::SameLine();
        if (ImGui::ArrowButton("##down", ImGuiDir_Down) && i + 1 < (int)cfg.sort.size()) {
            std::swap(cfg.sort[i], cfg.sort[i + 1]);
            sort_changed = true;
        }
        ImGui::PopID();
    }
    if (sort_changed) rebuild_order();

    ImGui::SeparatorText("Files");
    ImGui::Text("Settings are kept in %s", config_dir().c_str());
    ImGui::EndChild();
}

}  // namespace msv
