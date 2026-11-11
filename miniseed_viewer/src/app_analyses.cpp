// The Analyses tab: new traces (or functions) made from pairs of traces chosen on the Main tab.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "app.h"
#include "implot.h"

namespace msv {

namespace {

std::string shift_text(double s) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+.6g s", s);
    return buf;
}

}  // namespace

int App::find_trace(const std::string& id) const {
    for (size_t i = 0; i < traces_.size(); i++)
        if (traces_[i].id == id) return (int)i;
    return -1;
}

void App::remove_derived(int ti) {
    if (ti < 0 || ti >= (int)traces_.size() || !traces_[ti].derived) return;
    derived_calib_.erase(traces_[ti].id);
    traces_.erase(traces_.begin() + ti);
    topts_.erase(topts_.begin() + ti);
    calib_sel_.erase(calib_sel_.begin() + ti);
    last_used_.erase(last_used_.begin() + ti);
    clear_caches();  // everything cached is filed under trace numbers, which have just shifted
    rebuild_order();
}

// Prepares the two chosen traces (reading them from disk if needed). Returns false after setting the message.
const FilterSpec* App::analysis_filter() const {
    return an_.filter >= 0 && an_.filter < (int)filters_.size() ? &filters_[an_.filter] : nullptr;
}

// f1 and f2 are the filtered copies (empty when no filter is chosen, or when none is needed: the difference
// only needs them to find R). Filtering is done on the whole record, then the window is cut out, so the
// filter's start-up effects stay outside the window.
bool App::analysis_inputs(Pieces& p1, Pieces& p2, Pieces& f1, Pieces& f2, int& t1, int& t2) {
    an_.message_is_error = true;
    t1 = find_trace(an_.in1);
    t2 = find_trace(an_.in2);
    if (t1 < 0 || t2 < 0) {
        an_.message = "Choose two traces first.";
        return false;
    }
    std::vector<Trace*> need;
    for (int ti : {t1, t2})
        if (!traces_[ti].loaded && !traces_[ti].derived) need.push_back(&traces_[ti]);
    if (!need.empty()) {
        if (need.size() == 2 && need[0] == need[1]) need.pop_back();
        ::msv::load_traces(need, messages_);
    }
    const bool win_view = an_.kind >= 1 && an_.use_view && !an_.view_sig.empty();
    const bool win = win_view || (an_.kind == 0 && an_.use_sel && sel_state_ == 2);
    // In "seconds since start" mode the window is converted to clock time with trace 1's start, and the same
    // clock window is used for both traces, so the two always line up by clock time.
    const double woff = off(t1);
    const double inf = std::numeric_limits<double>::infinity();
    auto make = [&](int ti) {
        auto sh = an_.shifts.find(traces_[ti].id);
        double shift = sh == an_.shifts.end() ? 0.0 : sh->second;
        return prepare_pieces(traces_[ti], calib_for(ti).scale, flips_.count(traces_[ti].id) > 0, shift, -inf, inf);
    };
    p1 = make(t1);
    p2 = make(t2);
    std::string err;
    if (!match_rates(p1, p2, an_.rate_note, err)) {
        an_.message = err;
        return false;
    }
    f1.clear();
    f2.clear();
    const FilterSpec* filt = analysis_filter();
    if (filt && (an_.kind >= 1 || an_.balance)) {
        const double rate = 1.0 / p1[0].dt;
        std::string why = check_filter(*filt, rate);
        if (!why.empty()) {
            an_.message = "The filter cannot be used at this sample rate (" + std::to_string(rate) + "/s): " + why;
            return false;
        }
        f1 = p1;
        f2 = p2;
        filter_pieces(f1, *filt);
        filter_pieces(f2, *filt);
    }
    if (win) {
        const double w0 = win_view ? an_.view[0] : sel_[0], w1 = win_view ? an_.view[1] : sel_[1];
        for (Pieces* ps : {&p1, &p2, &f1, &f2})
            if (!ps->empty()) *ps = cut_pieces(*ps, woff + w0, woff + w1);
    }
    return true;
}

void App::run_diff() {
    Pieces p1, p2, f1, f2;
    int t1, t2;
    if (!analysis_inputs(p1, p2, f1, f2, t1, t2)) return;
    std::vector<Segment> segs;
    std::string err;
    const FilterSpec* filt = an_.balance ? analysis_filter() : nullptr;
    if (!diff_pieces(p1, p2, an_.balance, &an_.balance_R, segs, err, filt ? &f1 : nullptr, filt ? &f2 : nullptr)) {
        an_.message = err;
        return;
    }
    Trace t;
    t.derived = true;
    t.net = "DIFF";
    t.sta = std::to_string(++an_.serial);
    t.id = t.net + t.sta;
    t.note = traces_[t1].id + " - " + traces_[t2].id;
    for (int ti : {t1, t2}) {
        auto sh = an_.shifts.find(traces_[ti].id);
        if (sh != an_.shifts.end() && sh->second != 0) t.note += "  [" + traces_[ti].id + " shifted " + shift_text(sh->second) + "]";
    }
    if (an_.use_sel && sel_state_ == 2) t.note += "  [Main window only]";
    if (an_.balance) {
        char b[64];
        std::snprintf(b, sizeof b, "  [R = %.6g: sqrt(R)*tr1 - tr2/sqrt(R)]", an_.balance_R);
        t.note += b;
        if (filt) t.note += "  [R found after filter " + filt->name + "]";
    }
    if (!an_.rate_note.empty()) t.note += "  [" + an_.rate_note + "]";
    t.file_order = (1 << 20) + an_.serial;
    t.segs = std::move(segs);
    t.v_min = std::numeric_limits<double>::infinity();
    t.v_max = -t.v_min;
    bool first = true;
    for (Segment& s : t.segs) {
        s.build_summary();
        double mn, mx;
        s.minmax(0, s.size(), mn, mx);
        t.v_min = std::min(t.v_min, mn);
        t.v_max = std::max(t.v_max, mx);
        t.nsamples += s.size();
        if (first) {
            t.t_min = s.t0;
            t.rate = 1.0 / s.dt;
            first = false;
        }
        t.t_max = s.t_end();
    }
    t.loaded = true;
    const std::string u1 = calib_for(t1).units, u2 = calib_for(t2).units;
    Calib c;
    c.units = u1 == u2 ? u1 : "mixed units";
    derived_calib_[t.id] = c;
    const std::string id = t.id, note = t.note;
    traces_.push_back(std::move(t));
    topts_.resize(traces_.size());
    calib_sel_.resize(traces_.size(), 0);
    last_used_.resize(traces_.size(), 0);
    rebuild_order();
    const int n = std::max(1, cfg.n_per_screen);
    for (size_t k = 0; k < order_.size(); k++)
        if (traces_[order_[k]].id == id) page_ = (int)k / n;  // show the page that holds the new trace
    an_.message_is_error = false;
    an_.message = "Made trace " + id + " = " + note + ". It is also on the Main tab and can be plotted in the other tabs.";
    if (u1 != u2) an_.message += "  Warning: the two traces have different units (" + u1 + ", " + u2 + ").";
}

void App::run_xcorr() {
    Pieces p1, p2, f1, f2;
    int t1, t2;
    if (!analysis_inputs(p1, p2, f1, f2, t1, t2)) return;
    XcorrResult R;
    std::string err;
    const FilterSpec* filt = analysis_filter();
    if (!cross_correlate(filt ? f1 : p1, filt ? f2 : p2, an_.lag[0], an_.lag[1], an_.normalize, R, err)) {
        an_.message = err;
        return;
    }
    an_.xc = std::move(R);
    an_.have_xc = true;
    an_.xc_refit = true;
    an_.xc_normalized = an_.normalize;
    an_.xc_title = traces_[t1].id + " x " + traces_[t2].id;
    an_.message_is_error = false;
    an_.message = "Cross-correlated " + std::to_string(an_.xc.n_samples) + " samples (" +
                  std::to_string(an_.xc.n_samples * an_.xc.dt) + " s)." +
                  (an_.rate_note.empty() ? "" : "  Note: " + an_.rate_note + ".") +
                  (filt ? "  Filter: " + filt->name + "." : "");
    if (an_.use_view && !an_.view_sig.empty()) an_.message += "  Interval: " + fmt_disp(an_.view[0]) + " to " + fmt_disp(an_.view[1]) + ".";
}

void App::run_correlogram() {
    Pieces p1, p2, f1, f2;
    int t1, t2;
    if (!analysis_inputs(p1, p2, f1, f2, t1, t2)) return;
    CorrelogramResult C;
    std::string err;
    const FilterSpec* filt = analysis_filter();
    if (!correlogram(filt ? f1 : p1, filt ? f2 : p2, an_.cg_win, an_.cg_overlap, an_.lag[0], an_.lag[1], C, err)) {
        an_.message = err;
        return;
    }
    auto D = std::make_shared<SgData>();
    D->w = C.ncols;
    D->h = C.nlags;
    const double woff = off(t1);  // the plots' time axis: clock time, or seconds since trace 1 started
    D->t0 = C.t_first + 0.5 * C.win - 0.5 * C.hop - woff;
    D->t1 = C.t_first + (C.ncols - 1) * C.hop + 0.5 * C.win + 0.5 * C.hop - woff;
    D->f0 = C.lag0 - 0.5 * C.dt;
    D->f1 = C.lag0 + (C.nlags - 1) * C.dt + 0.5 * C.dt;
    D->db = std::move(C.c);
    {   // typical size of the coefficients: 99.5th percentile of |value| (a sample of at most ~50000 values)
        std::vector<float> a;
        const size_t stride = std::max<size_t>(1, D->db.size() / 50000);
        for (size_t i = 0; i < D->db.size(); i += stride)
            if (!std::isnan(D->db[i])) a.push_back(std::fabs(D->db[i]));
        if (!a.empty()) {
            size_t k = std::min(a.size() - 1, (size_t)(0.995 * a.size()));
            std::nth_element(a.begin(), a.begin() + k, a.end());
            D->p_hi = std::min(1.0, (double)a[k]);
        } else {
            D->p_hi = 1.0;
        }
    }
    int n_ok = 0;
    for (int k = 0; k < D->w; k++)
        if (!std::isnan(D->db[(size_t)(D->h / 2) * D->w + k])) n_ok++;
    an_.cg = D;
    an_.have_cg = true;
    an_.cg_title = traces_[t1].id + " x " + traces_[t2].id;
    an_.message_is_error = false;
    an_.message = "Correlogram: " + std::to_string(D->w) + " windows of " + std::to_string(C.win) + " s (" +
                  std::to_string(n_ok) + " with data), " + std::to_string(D->h) + " lags." +
                  (an_.rate_note.empty() ? "" : "  Note: " + an_.rate_note + ".") +
                  (filt ? "  Filter: " + filt->name + "." : "");
    if (an_.use_view && !an_.view_sig.empty()) an_.message += "  Interval: " + fmt_disp(an_.view[0]) + " to " + fmt_disp(an_.view[1]) + ".";
}

void App::draw_analyses() {
    ImGui::BeginChild("##analyses", ImVec2(0, 0));
    draw_analyses_body();
    ImGui::EndChild();
}

// The two chosen traces, stacked, as they will go into the analysis: calibrated, polarity-corrected and shifted.
// Zooming and panning move the time axis only; each trace's own vertical scale fits what is on screen.
void App::draw_analysis_traces(float height) {
    const int ids[2] = {find_trace(an_.in1), find_trace(an_.in2)};
    if (ids[0] < 0 || ids[1] < 0) return;
    for (int ti : ids) {
        if (!traces_[ti].loaded && !traces_[ti].derived) {
            std::vector<Trace*> one{&traces_[ti]};
            ::msv::load_traces(one, messages_);
        }
        last_used_[ti] = frame_no_;  // keep it in memory
    }
    double shift[2], lo = 1e300, hi = -1e300;
    for (int k = 0; k < 2; k++) {
        auto sh = an_.shifts.find(traces_[ids[k]].id);
        shift[k] = sh == an_.shifts.end() ? 0.0 : sh->second;
        lo = std::min(lo, traces_[ids[k]].t_min + shift[k] - off(ids[0]));
        hi = std::max(hi, traces_[ids[k]].t_max + shift[k] - off(ids[0]));
    }
    // Fit the time axis to the data when the choice of traces (or the time axis type) changes.
    const std::string sig = an_.in1 + "|" + an_.in2 + "|" + (cfg.relative_time ? "r" : "c");
    if (sig != an_.view_sig || !(an_.view[1] > an_.view[0])) {
        an_.view_sig = sig;
        an_.view[0] = lo;
        an_.view[1] = hi > lo ? hi : lo + 1;
    }
    if (ImGui::SmallButton("Show all")) {
        an_.view[0] = lo;
        an_.view[1] = hi > lo ? hi : lo + 1;
    }
    if (sel_state_ == 2) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Zoom to the window selected on Main")) {
            an_.view[0] = sel_[0];
            an_.view[1] = sel_[1];
        }
    }
    ImGui::SameLine();
    if (an_.pick_state == 1)
        ImGui::Text("New t1 = %s. Click again to set t2 and zoom.", fmt_disp(an_.pick_first).c_str());
    else
        ImGui::TextDisabled("Click twice on a plot to zoom to a new interval (the wheel zooms too, right-drag pans). Traces are shown filtered and shifted as they will be analysed.");

    const bool cg = an_.kind == 2 && an_.have_cg && an_.cg && an_.cg->w > 0;
    const float bar_w = ImGui::GetFontSize() * 5.5f;
    if (cg) height = std::max(150.0f, ImGui::GetContentRegionAvail().y - 4);  // what is left below the buttons above
    if (cg) ImGui::BeginChild("##anplots", ImVec2(ImGui::GetContentRegionAvail().x - bar_w, height), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    const float label_w = ImGui::GetFontSize() * 11.0f;
    const float row_h = cg ? std::max(50.0f, std::floor(0.2f * (height - axis_label_height() - ImGui::GetTextLineHeightWithSpacing())))
                           : std::max(60.0f, std::floor((height - axis_label_height() - ImGui::GetTextLineHeightWithSpacing()) / 2));
    const int nbins = std::max(100, std::min(4000, (int)(ImGui::GetContentRegionAvail().x - label_w - 60)));
    bool clicked = false;
    double click_x = 0;
    auto mark_and_click = [&]() {  // the first of the two zoom clicks is marked; a click without dragging picks a time
        if (an_.pick_state == 1) {
            ImPlotSpec mark;
            mark.LineColor = mark_color();
            mark.LineWeight = 1.5f;
            mark.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##pick", &an_.pick_first, 1, mark);
        }
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            ImVec2 d = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
            if (d.x == 0 && d.y == 0) {
                clicked = true;
                click_x = ImPlot::GetPlotMousePos().x;
            }
        }
    };
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
    bool aligned = ImPlot::BeginAlignedPlots("##align_an");
    for (int k = 0; k < 2; k++) {
        const int ti = ids[k];
        const Trace& tr = traces_[ti];
        const bool flip = flips_.count(tr.id) > 0;
        double scale = calib_for(ti).scale * (flip ? -1.0 : 1.0);
        // With a filter chosen, plot the filtered trace (calibration, polarity and filter are then already in it).
        const std::vector<Segment>* src = &tr.segs;
        const FilterSpec* filt = analysis_filter();
        std::string fwarn;
        if (filt) {
            std::string why = check_filter(*filt, tr.rate);
            if (!why.empty()) {
                fwarn = "filter not applied: " + why;
            } else {
                char kb[320];
                std::snprintf(kb, sizeof kb, "%s|%.9g|%d|%.6f|%lld|%d,%d,%.9g,%.9g,%d", tr.id.c_str(), calib_for(ti).scale, (int)flip,
                              tr.t_min, (long long)tr.nsamples, (int)filt->type, filt->order, filt->f1, filt->f2, (int)filt->zero_phase);
                if (an_.pc_key[k] != kb) {
                    const double inf = std::numeric_limits<double>::infinity();
                    Pieces pc = prepare_pieces(tr, calib_for(ti).scale, flip, 0, -inf, inf);
                    filter_pieces(pc, *filt);
                    an_.pc_segs[k].clear();
                    for (const Piece& p : pc) {
                        Segment s;
                        s.t0 = p.t0;
                        s.dt = p.dt;
                        s.set(p.v);
                        s.build_summary();
                        an_.pc_segs[k].push_back(std::move(s));
                    }
                    an_.pc_key[k] = kb;
                }
                src = &an_.pc_segs[k];
                scale = 1.0;
            }
        }
        std::vector<Series> data;
        double vlo = 1e300, vhi = -1e300;
        const double shift_k = shift[k], o = off(ids[0]);  // both traces on one time axis (trace 1's start in relative mode)
        for (const Segment& seg : *src) {
            Series ser;
            seg.envelope(o + an_.view[0] - shift_k, o + an_.view[1] - shift_k, nbins, ser);
            if (ser.x.empty()) continue;
            for (size_t i = 0; i < ser.x.size(); i++) {
                ser.y[i] *= scale;
                vlo = std::min(vlo, ser.y[i]);
                vhi = std::max(vhi, ser.y[i]);
                ser.x[i] += shift_k - o;
            }
            data.push_back(std::move(ser));
        }
        if (vlo > vhi) vlo = 0, vhi = 1;
        double span = vhi - vlo, big = std::max(std::fabs(vlo), std::fabs(vhi));
        if (span <= 1e-12 * big || big == 0) {
            double m = big > 0 ? big : 1.0;
            vlo -= m;
            vhi += m;
        } else {
            vlo -= 0.05 * span;
            vhi += 0.05 * span;
        }
        ImGui::PushID(k);
        ImGui::BeginChild("lbl", ImVec2(label_w, row_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
        ImGui::TextColored(trace_color(ti), "%s", tr.id.c_str());
        ImGui::TextUnformatted(k == 0 ? "trace 1" : "trace 2");
        if (shift_k != 0) ImGui::TextUnformatted(shift_text(shift_k).c_str());
        if (filt && fwarn.empty()) ImGui::TextUnformatted(("filter: " + filt->name).c_str());
        if (!fwarn.empty()) ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.35f, 1), "%s", fwarn.c_str());
        const std::string& u = calib_for(ti).units;
        if (u != "counts") ImGui::TextUnformatted(u.c_str());
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        const bool last = k == 1 && !cg;
        if (ImPlot::BeginPlot("##p", ImVec2(-1, last ? row_h + axis_label_height() : row_h),
                              ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoFrame)) {
            ImPlot::SetupAxis(ImAxis_X1, nullptr,
                              ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight |
                                  (last ? 0 : ImPlotAxisFlags_NoTickLabels));
            ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoHighlight | ImPlotAxisFlags_Lock);
            ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
            ImPlot::SetupAxisLinks(ImAxis_X1, &an_.view[0], &an_.view[1]);
            ImPlot::SetupAxisLimits(ImAxis_Y1, vlo, vhi, ImPlotCond_Always);
            ImPlotSpec spec;
            spec.LineColor = trace_color(ti);
            spec.Flags = ImPlotItemFlags_NoLegend;
            for (const Series& s : data) ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);
            mark_and_click();
            ImPlot::EndPlot();
        }
        ImGui::PopID();
    }
    if (cg) {  // the correlogram: time along the same axis as the traces, lag up the side
        SgData& D = *an_.cg;
        // Colour range: symmetric about 0 (so 0 stays white). Automatic = the 99.5th percentile of |value|.
        const double lim = std::max(1e-6, an_.cg_auto ? D.p_hi : an_.cg_lim);
        char pbuf[64];
        std::snprintf(pbuf, sizeof pbuf, "cg|%.9g", lim);
        const std::string psig = pbuf;
        if (!D.tex || D.tex_sig != psig) make_picture(D, ImPlotColormap_RdBu, -lim, lim, psig);
        const float h = std::max(80.0f, ImGui::GetContentRegionAvail().y - 2);
        ImGui::BeginChild("lblcg", ImVec2(label_w, h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
        ImGui::TextUnformatted("correlogram");
        ImGui::TextUnformatted(an_.cg_title.c_str());
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        if (ImPlot::BeginPlot("##cgplot", ImVec2(-1, h),
                              ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame)) {
            ImPlot::SetupAxis(ImAxis_X1, nullptr, ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight);
            ImPlot::SetupAxis(ImAxis_Y1, "lag (s)", ImPlotAxisFlags_NoHighlight | ImPlotAxisFlags_Lock);  // only time zooms
            ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
            ImPlot::SetupAxisLinks(ImAxis_X1, &an_.view[0], &an_.view[1]);
            ImPlot::SetupAxisLimits(ImAxis_Y1, D.f0, D.f1, ImPlotCond_Always);
            {   // grey backdrop, so stretches without a result read as "no data"
                ImPlotRect lim = ImPlot::GetPlotLimits();
                ImVec2 p0 = ImPlot::PlotToPixels(lim.X.Min, lim.Y.Max), p1 = ImPlot::PlotToPixels(lim.X.Max, lim.Y.Min);
                ImPlot::PushPlotClipRect();
                ImPlot::GetPlotDrawList()->AddRectFilled(p0, p1, IM_COL32(80, 80, 88, 110));
                ImPlot::PopPlotClipRect();
            }
            if (D.tex)
                ImPlot::PlotImage("##img", (ImTextureID)(intptr_t)D.tex, ImPlotPoint(D.t0, D.f0), ImPlotPoint(D.t1, D.f1),
                                  ImVec2(0, 1), ImVec2(1, 0));
            mark_and_click();
            ImPlot::EndPlot();
        }
    }
    if (aligned) ImPlot::EndAlignedPlots();
    ImGui::PopStyleVar();
    if (clicked) {  // two clicks zoom to the interval between them (as in the Magnify tab)
        if (an_.pick_state == 0) {
            an_.pick_first = click_x;
            an_.pick_state = 1;
        } else {
            an_.pick_state = 0;
            double a = std::min(an_.pick_first, click_x), b = std::max(an_.pick_first, click_x);
            if (b - a >= 1e-6) {
                an_.view[0] = a;
                an_.view[1] = b;
            }
        }
    }
    if (cg) {
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        ImGui::BeginChild("##anbar", ImVec2(bar_w, height), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
        ImGui::Dummy(ImVec2(0, 4));
        const double lim = std::max(1e-6, an_.cg_auto ? an_.cg->p_hi : an_.cg_lim);
        ImPlot::ColormapScale("correlation coefficient", -lim, lim,
                              ImVec2(bar_w - 8, std::max(60.0f, height - axis_label_height() - 8)), "%.3g", 0, ImPlotColormap_RdBu);
        ImGui::EndChild();
    }
}

void App::draw_analyses_body() {
    const float fh = ImGui::GetFontSize();
    std::vector<int> cands;  // traces that are ticked on the Main tab (on any page)
    for (int ti : order_)
        if (traces_[ti].on) cands.push_back(ti);
    if (cands.empty()) {
        ImGui::Spacing();
        ImGui::TextWrapped("No traces are ticked on the Main tab. Tick two or more traces there, then come back here.");
        return;
    }
    auto listed = [&](const std::string& id) {
        int ti = find_trace(id);
        return ti >= 0 && traces_[ti].on;
    };
    if (!listed(an_.in1)) an_.in1 = traces_[cands[0]].id;
    if (!listed(an_.in2)) an_.in2 = traces_[cands[cands.size() > 1 ? 1 : 0]].id;

    ImGui::Spacing();
    combo_filter("Filter", an_.filter);
    ImGui::SameLine();
    if (ImGui::Button("Edit filters...")) {
        if (an_.filter >= 0) edit_filter_ = an_.filter;
        ImGui::OpenPopup("Filters");
    }
    filter_editor();
    ImGui::SameLine();
    ImGui::TextDisabled("Used for the trace plots below, cross-correlation, correlogram, and (for finding R) the difference.");
    ImGui::TextUnformatted("Operation:");
    ImGui::SameLine();
    ImGui::RadioButton("difference (trace 1 - trace 2)", &an_.kind, 0);
    ImGui::SameLine();
    ImGui::RadioButton("cross-correlation", &an_.kind, 1);
    ImGui::SameLine();
    ImGui::RadioButton("correlogram", &an_.kind, 2);

    auto trace_combo = [&](const char* label, std::string& id) {
        ImGui::SetNextItemWidth(fh * 16);
        if (ImGui::BeginCombo(label, id.c_str())) {
            for (int ti : cands)
                if (ImGui::Selectable((traces_[ti].id + "##c" + std::to_string(ti)).c_str(), traces_[ti].id == id)) id = traces_[ti].id;
            ImGui::EndCombo();
        }
    };
    trace_combo("trace 1", an_.in1);
    ImGui::SameLine(0, 20);
    trace_combo("trace 2", an_.in2);
    ImGui::SameLine();
    if (ImGui::SmallButton("swap")) std::swap(an_.in1, an_.in2);

    int nshift = 0;
    for (const auto& kv : an_.shifts)
        if (kv.second != 0) nshift++;
    const std::string hdr = "Time shifts (optional): " + std::to_string(nshift) + " set###shifts";
    if (ImGui::CollapsingHeader(hdr.c_str())) {
        ImGui::TextDisabled(
            "A positive shift moves a trace later in time (its data are treated as if recorded that many seconds later). "
            "Applied to every analysis below.");
        for (int ti : cands) {
            const std::string& id = traces_[ti].id;
            double v = 0;
            auto it = an_.shifts.find(id);
            if (it != an_.shifts.end()) v = it->second;
            ImGui::PushID(ti);
            ImGui::SetNextItemWidth(fh * 9);
            if (ImGui::InputDouble("##shift", &v, 0, 0, "%.6g")) {
                if (v == 0) an_.shifts.erase(id);
                else an_.shifts[id] = v;
            }
            ImGui::SameLine();
            ImGui::TextUnformatted(("s    " + id).c_str());
            ImGui::PopID();
        }
        if (ImGui::SmallButton("Clear all shifts")) an_.shifts.clear();
    }

    if (an_.kind == 0) {
        ImGui::BeginDisabled(sel_state_ != 2);
        ImGui::Checkbox("Use only the time window selected on the Main tab", &an_.use_sel);
        ImGui::EndDisabled();
        if (sel_state_ != 2) {
            ImGui::SameLine();
            ImGui::TextDisabled("(no window selected: whole traces are used)");
        }
    } else {
        ImGui::Checkbox("Use only the time interval shown in the plots below (zoom there to choose it)", &an_.use_view);
    }
    if (an_.kind == 2) {
        ImGui::SetNextItemWidth(fh * 9);
        ImGui::InputDouble("window length (s)", &an_.cg_win, 0, 0, "%.6g");
        ImGui::SameLine(0, 20);
        ImGui::SetNextItemWidth(fh * 6);
        ImGui::InputInt("overlap (%)", &an_.cg_overlap, 0, 0);
        an_.cg_overlap = std::max(0, std::min(95, an_.cg_overlap));
        ImGui::SameLine(0, 20);
    }

    if (an_.kind == 0)
        ImGui::Checkbox("Find the scale factor R that makes the difference smallest: output = sqrt(R)*trace 1 - trace 2/sqrt(R)",
                        &an_.balance);
    if (an_.kind == 0 && an_.balance) {
        ImGui::Indent(fh * 2);
        ImGui::TextDisabled("(if a filter is chosen at the top, it is applied to both traces to find R only; the output is never filtered)");
        ImGui::Unindent(fh * 2);
    }
    if (an_.kind >= 1) {
        ImGui::SetNextItemWidth(fh * 9);
        ImGui::InputDouble("smallest lag (s)", &an_.lag[0], 0, 0, "%.6g");
        ImGui::SameLine(0, 20);
        ImGui::SetNextItemWidth(fh * 9);
        ImGui::InputDouble("largest lag (s)", &an_.lag[1], 0, 0, "%.6g");
        if (an_.kind == 1) {
            ImGui::SameLine(0, 20);
            ImGui::Checkbox("normalize (1 = identical)", &an_.normalize);
        }
    }
    if (an_.kind == 2) {
        ImGui::Checkbox("autoscale color bar", &an_.cg_auto);
        if (!an_.cg_auto) {
            ImGui::SameLine(0, 20);
            ImGui::SetNextItemWidth(fh * 7);
            ImGui::InputDouble("color range: -limit to +limit", &an_.cg_lim, 0, 0, "%.3g");
            an_.cg_lim = std::max(an_.cg_lim, 1e-6);
        }
    }
    if (ImGui::Button(an_.kind == 0 ? "Make difference trace" : "Calculate")) {
        if (an_.kind == 0) run_diff();
        else if (an_.kind == 1) run_xcorr();
        else run_correlogram();
    }
    if (!an_.message.empty()) {
        ImGui::SameLine(0, 20);
        ImGui::PushStyleColor(ImGuiCol_Text, an_.message_is_error ? ImVec4(1.0f, 0.55f, 0.35f, 1) : ImGui::GetStyleColorVec4(ImGuiCol_Text));
        ImGui::PushTextWrapPos(0);
        ImGui::TextUnformatted(an_.message.c_str());
        ImGui::PopTextWrapPos();
        ImGui::PopStyleColor();
    }

    // Traces made here, with a way to delete them.
    int n_derived = 0;
    for (const Trace& t : traces_) n_derived += t.derived;
    if (n_derived > 0) {
        const std::string h = "Traces made by analyses (" + std::to_string(n_derived) + ")###derived";
        if (ImGui::CollapsingHeader(h.c_str())) {
            int del = -1;
            for (int ti = 0; ti < (int)traces_.size(); ti++) {
                if (!traces_[ti].derived) continue;
                ImGui::PushID(ti);
                if (ImGui::SmallButton("delete")) del = ti;
                ImGui::SameLine();
                ImGui::Text("%s = %s", traces_[ti].id.c_str(), traces_[ti].note.c_str());
                ImGui::PopID();
            }
            if (del >= 0) remove_derived(del);
        }
    }

    {
        const bool xc = an_.have_xc && an_.kind == 1;
        const float avail = ImGui::GetContentRegionAvail().y;
        ImGui::Separator();
        draw_analysis_traces(xc ? std::max(220.0f, avail * 0.45f) : std::max(220.0f, avail - 4));
    }

    if (an_.have_xc && an_.kind == 1) {
        const XcorrResult& R = an_.xc;
        ImGui::Separator();
        ImGui::Text("%s:  largest |value| %.4g at lag %+.4f s.  Positive lag: trace 2 is later than trace 1.",
                    an_.xc_title.c_str(), R.peak_value, R.peak_lag);
        ImGui::SameLine(0, 20);
        if (ImGui::SmallButton("Fit to screen")) an_.xc_refit = true;
        if (an_.xc_refit) {
            ImPlot::SetNextAxesToFit();
            an_.xc_refit = false;
        }
        if (ImPlot::BeginPlot("##xc", ImVec2(-1, std::max(120.0f, ImGui::GetContentRegionAvail().y - 4)),
                              ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect)) {
            ImPlot::SetupAxes("lag (s)", an_.xc_normalized ? "correlation" : "sum of products");
            ImPlotSpec spec;
            spec.LineColor = trace_color(0);
            spec.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotLine("##xcl", R.lag.data(), R.c.data(), (int)R.lag.size(), spec);
            ImPlotSpec mark;
            mark.LineColor = mark_color();
            mark.LineWeight = 1.5f;
            mark.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##peak", &R.peak_lag, 1, mark);
            ImPlot::EndPlot();
        }
    }
}

}  // namespace msv
