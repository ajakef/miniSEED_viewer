// The spectrogram tab: power versus time and frequency for each selected trace.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdint>

#include <GLFW/glfw3.h>  // for the OpenGL calls that make the pictures

#include "app.h"
#include "imgui.h"
#include "implot.h"
#include "parallel.h"
#include "pocketfft_hdronly.h"

#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif

namespace msv {

namespace {

const ImPlotColormap kCmaps[] = {ImPlotColormap_Viridis, ImPlotColormap_Plasma, ImPlotColormap_Hot,
                                 ImPlotColormap_Jet,     ImPlotColormap_Spectral, ImPlotColormap_Cool};
const char* kCmapNames = "Viridis\0Plasma\0Hot\0Jet\0Spectral\0Cool\0";

// Spectrogram of the processed data between absolute times v1 and v2.
// Each column of the result averages the power of all windows centred in it.
std::shared_ptr<SgData> compute_sg(const Processed& P, double v1, double v2, double off, int nwin, int overlap_pct,
                                   bool log_f) {
    auto R = std::make_shared<SgData>();
    nwin = std::max(32, std::min(16384, nwin));
    struct Piece {
        const Segment* s;
        int64_t i0, i1;
    };
    std::vector<Piece> pieces;
    double dt = 0;
    for (const Segment& s : P.segs) {
        int64_t i0, i1;
        s.index_range(v1, v2, i0, i1);
        if (i1 - i0 < nwin) continue;
        if (dt == 0) dt = s.dt;
        else if (std::fabs(s.dt - dt) > 1e-6 * dt) continue;
        pieces.push_back({&s, i0, i1});
    }
    if (pieces.empty()) return R;

    const double fs = 1.0 / dt;
    const int bins = nwin / 2 + 1;
    int ncols = std::max(16, std::min(2048, 4000000 / bins));
    int64_t hop = std::max<int64_t>(1, (int64_t)nwin * (100 - std::max(0, std::min(95, overlap_pct))) / 100);
    // Do not compute many more windows than the picture can show.
    int64_t col_samples = (int64_t)((v2 - v1) / dt / ncols);
    hop = std::max<int64_t>(hop, col_samples / 8);

    int64_t nframes = 0;
    double tfirst = 1e300, tlast = -1e300;
    for (const Piece& pc : pieces) {
        int64_t nf = (pc.i1 - pc.i0 - nwin) / hop + 1;
        nframes += nf;
        tfirst = std::min(tfirst, pc.s->time_at(pc.i0 + nwin / 2));
        tlast = std::max(tlast, pc.s->time_at(pc.i0 + (nf - 1) * hop + nwin / 2));
    }
    ncols = (int)std::max<int64_t>(1, std::min<int64_t>(ncols, nframes));
    double ta = tfirst, tb = tlast;
    if (tb <= ta) {  // a single window
        ta -= 0.5 * hop * dt;
        tb += 0.5 * hop * dt;
    }
    const double col_dt = (tb - ta) / ncols;

    std::vector<double> win(nwin);
    double wsum2 = 0;
    for (int i = 0; i < nwin; i++) {
        win[i] = 0.5 - 0.5 * std::cos(2 * 3.14159265358979323846 * i / (nwin - 1));
        wsum2 += win[i] * win[i];
    }
    const double scale = 2.0 / (fs * wsum2);

    std::vector<double> acc((size_t)ncols * bins, 0.0);
    std::vector<int> cnt(ncols, 0);
    std::vector<double> buf(nwin);
    std::vector<std::complex<double>> X(bins);
    for (const Piece& pc : pieces) {
        int64_t nf = (pc.i1 - pc.i0 - nwin) / hop + 1;
        for (int64_t k = 0; k < nf; k++) {
            int64_t i = pc.i0 + k * hop;
            double mean = 0;
            for (int j = 0; j < nwin; j++) {
                buf[j] = pc.s->at(i + j);
                mean += buf[j];
            }
            mean /= nwin;
            for (int j = 0; j < nwin; j++) buf[j] = (buf[j] - mean) * win[j];
            pocketfft::r2c(pocketfft::shape_t{(size_t)nwin}, pocketfft::stride_t{(ptrdiff_t)sizeof(double)},
                           pocketfft::stride_t{(ptrdiff_t)sizeof(std::complex<double>)}, 0, true, buf.data(), X.data(),
                           1.0);
            double tc = pc.s->time_at(i + nwin / 2);
            int c = std::max(0, std::min(ncols - 1, (int)((tc - ta) / col_dt)));
            double* a = &acc[(size_t)c * bins];
            for (int b = 0; b < bins; b++) a[b] += std::norm(X[b]) * scale * ((b == 0 || (nwin % 2 == 0 && b == bins - 1)) ? 0.5 : 1.0);
            cnt[c]++;
        }
    }

    // Frequency rows: evenly spaced on a linear axis, evenly spaced in log(f) on a log axis.
    const double df = fs / nwin, nyq = fs / 2;
    int nrows;
    std::vector<double> edges;
    if (log_f) {
        nrows = 300;
        edges.resize(nrows + 1);
        for (int r = 0; r <= nrows; r++) edges[r] = df * std::pow(nyq / df, (double)r / nrows);
    } else {
        nrows = std::min(bins - 1, 512);
        edges.resize(nrows + 1);
        for (int r = 0; r <= nrows; r++) edges[r] = nyq * (double)r / nrows;
    }
    struct RowMap {
        int k0, k1;        // bins [k0, k1) are averaged; if empty, interpolate at frac
        double frac;
    };
    std::vector<RowMap> rm(nrows);
    for (int r = 0; r < nrows; r++) {
        int k0 = (int)std::ceil(edges[r] / df - 1e-9), k1 = (int)std::ceil(edges[r + 1] / df - 1e-9);
        if (log_f && k0 < 1) k0 = 1;
        k1 = std::min(k1, bins);
        double kc = 0.5 * (edges[r] + edges[r + 1]) / df;
        rm[r] = {k0, k1, kc};
    }
    R->w = ncols;
    R->h = nrows;
    R->db.assign((size_t)ncols * nrows, std::nanf(""));
    std::vector<float> sample;
    const size_t stride = std::max<size_t>(1, (size_t)ncols * nrows / 50000);
    for (int c = 0; c < ncols; c++) {
        if (cnt[c] == 0) continue;
        const double* a = &acc[(size_t)c * bins];
        const double inv = 1.0 / cnt[c];
        for (int r = 0; r < nrows; r++) {
            double p;
            if (rm[r].k1 > rm[r].k0) {
                double sum = 0;
                for (int k = rm[r].k0; k < rm[r].k1; k++) sum += a[k];
                p = sum / (rm[r].k1 - rm[r].k0);
            } else {
                int k = std::min(bins - 2, std::max(0, (int)std::floor(rm[r].frac)));
                double fr = std::max(0.0, std::min(1.0, rm[r].frac - k));
                p = (1 - fr) * a[k] + fr * a[k + 1];
            }
            float v = (float)(10.0 * std::log10(std::max(p * inv, 1e-30)));
            size_t idx = (size_t)r * ncols + c;
            R->db[idx] = v;
            if (idx % stride == 0) sample.push_back(v);
        }
    }
    if (sample.empty()) {
        R->w = R->h = 0;
        return R;
    }
    auto pct = [&](double q) {
        size_t k = std::min(sample.size() - 1, (size_t)(q * sample.size()));
        std::nth_element(sample.begin(), sample.begin() + k, sample.end());
        return (double)sample[k];
    };
    R->p_lo = pct(0.02);
    R->p_hi = pct(0.995);
    if (R->p_hi < -250) {  // no signal at all (a constant trace)
        R->w = R->h = 0;
        R->db.clear();
        return R;
    }
    R->t0 = ta - off;
    R->t1 = tb - off;
    R->f0 = edges.front();
    R->f1 = edges.back();
    return R;
}

// Makes (or remakes) the picture for a spectrogram with the given colours.
void upload_picture(SgData& d, ImPlotColormap cmap, double dbmin, double dbmax, const std::string& sig) {
    if (d.w <= 0 || d.h <= 0) return;
    uint32_t lut[256];
    for (int i = 0; i < 256; i++) lut[i] = ImGui::ColorConvertFloat4ToU32(ImPlot::SampleColormap(i / 255.0f, cmap));
    const double span = std::max(dbmax - dbmin, 1e-9);
    std::vector<uint32_t> px(d.db.size());
    for (size_t i = 0; i < px.size(); i++) {
        float v = d.db[i];
        if (std::isnan(v)) {
            px[i] = 0;  // transparent: no data here
            continue;
        }
        int k = (int)(255.0 * (v - dbmin) / span);
        px[i] = lut[std::max(0, std::min(255, k))];
    }
    GLuint t = d.tex;
    if (!t) glGenTextures(1, &t);
    glBindTexture(GL_TEXTURE_2D, t);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, d.w, d.h, 0, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
    d.tex = t;
    d.tex_sig = sig;
}

}  // namespace

SgData::~SgData() {
    if (tex) {
        GLuint t = tex;
        glDeleteTextures(1, &t);
    }
}

// ---------------------------------------------------------------------------

std::vector<std::shared_ptr<SgData>> App::specgrams_for(const std::vector<int>& rows, bool wf, const SgState& sg) {
    std::vector<std::shared_ptr<SgData>> out(rows.size());
    const bool stable = view_stable();
    std::vector<std::string> keys(rows.size()), stale_keys(rows.size());
    std::vector<int> todo;
    for (size_t r = 0; r < rows.size(); r++) {
        int ti = rows[r];
        char buf[160];
        std::snprintf(buf, sizeof buf, "|v%.6f|%.6f|%d|%d|%d", off(ti) + view_[0], off(ti) + view_[1], sg.nwin,
                      sg.overlap_pct, (int)sg.log_f);
        keys[r] = proc_key(ti, wf) + buf;
        stale_keys[r] = std::to_string(ti) + "|" + std::to_string((int)wf);
        auto it = sg_cache_.find(keys[r]);
        if (it != sg_cache_.end()) {
            out[r] = it->second;
            sg_stale_[stale_keys[r]] = it->second;
        } else if (stable) {
            todo.push_back((int)r);
        } else {
            auto st = sg_stale_.find(stale_keys[r]);
            if (st != sg_stale_.end()) out[r] = st->second;
        }
    }
    if (todo.empty()) return out;
    if (sg_cache_.size() > (size_t)(2 * std::max(1, cfg.n_per_screen) + 8)) {
        // Too many pictures kept: start afresh with the rows on screen.
        sg_cache_.clear();
        sg_stale_.clear();
        std::fill(out.begin(), out.end(), nullptr);
        todo.clear();
        for (size_t r = 0; r < rows.size(); r++) todo.push_back((int)r);
    }
    std::vector<std::shared_ptr<SgData>> res(todo.size());
    parallel_for((int)todo.size(), [&](int k) {
        int ti = rows[todo[k]];
        auto P = proc_for(ti, wf);
        if (!P) return;
        res[k] = compute_sg(*P, off(ti) + view_[0], off(ti) + view_[1], off(ti), sg.nwin, sg.overlap_pct, sg.log_f);
    });
    for (size_t k = 0; k < todo.size(); k++) {
        if (!res[k]) continue;
        sg_cache_[keys[todo[k]]] = res[k];
        sg_stale_[stale_keys[todo[k]]] = res[k];
        out[todo[k]] = res[k];
    }
    return out;
}

void App::specgram_toolbar(SgState& sg) {
    ImGui::PushID(&sg);
    const float fh = ImGui::GetFontSize();
    ImGui::SetNextItemWidth(fh * 7);
    ImGui::InputInt("window (samples, then Enter)", &sg.nwin, 0, 0, ImGuiInputTextFlags_EnterReturnsTrue);
    sg.nwin = std::max(32, std::min(16384, sg.nwin));
    ImGui::SameLine(0, 20);
    ImGui::SetNextItemWidth(fh * 5);
    ImGui::InputInt("overlap % (then Enter)", &sg.overlap_pct, 0, 0, ImGuiInputTextFlags_EnterReturnsTrue);
    sg.overlap_pct = std::max(0, std::min(95, sg.overlap_pct));
    ImGui::SameLine(0, 20);
    ImGui::Checkbox("Log frequency", &sg.log_f);
    ImGui::SameLine(0, 20);
    ImGui::SetNextItemWidth(fh * 8);
    ImGui::Combo("Colors", &sg.cmap, kCmapNames);
    ImGui::SameLine(0, 20);
    ImGui::Checkbox("Same filter as Magnify tab", &sg.same_filter);
    if (!sg.same_filter) {
        ImGui::SameLine();
        ImGui::Checkbox("Filtered", &sg.filtered);
    }

    ImGui::Checkbox("Automatic color range", &sg.auto_db);
    const float w = fh * 5;
    if (!sg.auto_db) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(w);
        ImGui::InputDouble("dB min", &sg.db[0], 0, 0, "%.0f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(w);
        ImGui::InputDouble("dB max", &sg.db[1], 0, 0, "%.0f");
    }
    ImGui::SameLine(0, 20);
    ImGui::TextUnformatted("Frequency limits (0 = automatic):");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(w);
    ImGui::InputDouble("f min", &sg.flim[0], 0, 0, "%.3g");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(w);
    ImGui::InputDouble("f max", &sg.flim[1], 0, 0, "%.3g");
    ImGui::PopID();
}

void App::draw_specgram(const char* id, float height, SgState& sg) {
    ImGui::BeginChild(id, ImVec2(0, height), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    if (!have_win_) {
        need_window_message();
        ImGui::EndChild();
        return;
    }
    specgram_toolbar(sg);
    std::vector<int> rows = page_traces(true);
    if (rows.empty()) {
        ImGui::TextUnformatted("All traces on this page are turned off (see the Main tab).");
        ImGui::EndChild();
        return;
    }
    const bool wf = sg.same_filter ? true : sg.filtered;
    ensure_processed(rows, wf);
    auto res = specgrams_for(rows, wf, sg);
    const ImPlotColormap cmap = kCmaps[std::max(0, std::min(5, sg.cmap))];

    // Colour range: automatic = from the lowest typical to the highest typical power of all rows.
    double dbmin = sg.db[0], dbmax = sg.db[1];
    if (sg.auto_db) {
        dbmin = 1e300;
        dbmax = -1e300;
        for (const auto& r : res)
            if (r && r->w > 0) {
                dbmin = std::min(dbmin, r->p_lo);
                dbmax = std::max(dbmax, r->p_hi);
            }
        if (dbmin > dbmax) {
            dbmin = -100;
            dbmax = 0;
        }
        if (dbmax - dbmin < 1) dbmax = dbmin + 1;
    }
    char sigbuf[96];
    std::snprintf(sigbuf, sizeof sigbuf, "%d|%.3f|%.3f", sg.cmap, dbmin, dbmax);
    for (const auto& r : res)
        if (r && r->w > 0 && r->tex_sig != sigbuf) upload_picture(*r, cmap, dbmin, dbmax, sigbuf);

    // Shared frequency axis: reset when the choice of traces or mode changes.
    {
        std::string sig;
        for (int ti : rows) sig += traces_[ti].id + ",";
        sig += std::to_string((int)wf) + std::to_string(sg.nwin) + std::to_string((int)sg.log_f);
        double fmin = 1e300, fmax = -1e300;
        for (const auto& r : res)
            if (r && r->w > 0) {
                fmin = std::min(fmin, r->f0);
                fmax = std::max(fmax, r->f1);
            }
        double want[2] = {fmin, fmax};
        bool have = fmin <= fmax;
        if (sg.flim[0] >= 0 && sg.flim[1] > sg.flim[0] && (sg.flim[0] > 0 || sg.flim[1] > 0)) {
            want[0] = sg.log_f ? std::max(sg.flim[0], 1e-6) : sg.flim[0];
            want[1] = sg.flim[1];
            have = true;
        }
        if (have && (want[0] != sg.last_desired[0] || want[1] != sg.last_desired[1] || sig != sg.last_sig)) {
            sg.yv[0] = sg.last_desired[0] = want[0];
            sg.yv[1] = sg.last_desired[1] = want[1];
            sg.last_sig = sig;
        }
    }

    const float fh = ImGui::GetFontSize();
    const float label_w = fh * 11.0f;
    const float bar_w = fh * 5.5f;
    const float avail_h = ImGui::GetContentRegionAvail().y;
    const float row_h = std::max(30.0f, std::floor((avail_h - axis_label_height()) / (int)rows.size()));
    double dlo, dhi;
    data_range(dlo, dhi);
    const double dspan = dhi - dlo;
    const ImPlotFlags pflags =
        ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

    ImGui::BeginChild("##sgplots", ImVec2(ImGui::GetContentRegionAvail().x - bar_w, avail_h), ImGuiChildFlags_None,
                      ImGuiWindowFlags_NoScrollbar);
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
    bool aligned = ImPlot::BeginAlignedPlots("##align_sg");
    for (size_t r = 0; r < rows.size(); r++) {
        int ti = rows[r];
        ImGui::PushID(ti);
        ImGui::BeginChild("lbl", ImVec2(label_w, row_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
        const char* note = nullptr;
        if (res[r] && res[r]->w == 0) note = "no signal in view";
        trace_label(ti, false, false, note);
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        const bool last = (r + 1 == rows.size());
        if (ImPlot::BeginPlot("##p", ImVec2(-1, last ? row_h + axis_label_height() : row_h), pflags)) {
            ImPlot::SetupAxis(ImAxis_X1, nullptr,
                              ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight |
                                  (last ? 0 : ImPlotAxisFlags_NoTickLabels));
            ImPlot::SetupAxis(ImAxis_Y1, "Hz", ImPlotAxisFlags_NoHighlight);
            ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
            ImPlot::SetupAxisScale(ImAxis_Y1, sg.log_f ? ImPlotScale_Log10 : ImPlotScale_Linear);
            ImPlot::SetupAxisLinks(ImAxis_X1, &view_[0], &view_[1]);
            ImPlot::SetupAxisLinks(ImAxis_Y1, &sg.yv[0], &sg.yv[1]);
            ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, dlo - dspan, dhi + dspan);
            if (res[r] && res[r]->w > 0 && res[r]->tex)
                ImPlot::PlotImage("##img", (ImTextureID)(intptr_t)res[r]->tex, ImPlotPoint(res[r]->t0, res[r]->f0),
                                  ImPlotPoint(res[r]->t1, res[r]->f1), ImVec2(0, 1), ImVec2(1, 0));
            ImPlot::EndPlot();
        }
        ImGui::PopID();
    }
    if (aligned) ImPlot::EndAlignedPlots();
    ImGui::PopStyleVar();
    ImGui::EndChild();
    ImGui::SameLine(0, 0);
    ImGui::BeginChild("##sgbar", ImVec2(bar_w, avail_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    ImGui::Dummy(ImVec2(0, 4));
    const std::string u = units_of(rows);
    ImPlot::ColormapScale((u == "mixed units" ? std::string("dB (mixed units)") : "dB re 1 " + u + "\xc2\xb2/Hz").c_str(),
                          dbmin, dbmax,
                          ImVec2(bar_w - 8, std::max(60.0f, avail_h - axis_label_height() - 8)), "%.0f", 0, cmap);
    ImGui::EndChild();
    ImGui::EndChild();
}

}  // namespace msv
