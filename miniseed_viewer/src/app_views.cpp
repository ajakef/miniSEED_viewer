// Magnify, spectra and split tabs.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "app.h"
#include "imgui.h"
#include "implot.h"

namespace msv {

namespace {
const char* kFilterTypeNames = "low-pass\0high-pass\0band-pass\0band-stop\0";
}

// ---------------------------------------------------------------------------
// Filter window

void App::filter_editor() {
    if (!ImGui::BeginPopupModal("Filters", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;
    const float fh = ImGui::GetFontSize();
    edit_filter_ = std::max(0, std::min(edit_filter_, (int)filters_.size() - 1));
    if (!filters_.empty() && draft_src_ != edit_filter_) {
        draft_ = filters_[edit_filter_];
        draft_src_ = edit_filter_;
    }
    ImGui::BeginChild("flist", ImVec2(fh * 14, fh * 12), ImGuiChildFlags_Borders);
    for (int i = 0; i < (int)filters_.size(); i++)
        if (ImGui::Selectable((filters_[i].name + "##fl" + std::to_string(i)).c_str(), edit_filter_ == i))
            edit_filter_ = i;
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginGroup();
    char name[128];
    std::snprintf(name, sizeof name, "%s", draft_.name.c_str());
    ImGui::SetNextItemWidth(fh * 14);
    if (ImGui::InputText("name", name, sizeof name)) draft_.name = name;
    int type = (int)draft_.type;
    ImGui::SetNextItemWidth(fh * 8);
    if (ImGui::Combo("type", &type, kFilterTypeNames)) draft_.type = (FilterType)type;
    ImGui::SetNextItemWidth(fh * 6);
    ImGui::InputInt("order", &draft_.order);
    draft_.order = std::max(1, std::min(10, draft_.order));
    const bool band = draft_.type == FilterType::BandPass || draft_.type == FilterType::BandStop;
    ImGui::SetNextItemWidth(fh * 8);
    ImGui::InputDouble(band ? "lower corner (Hz)" : "corner (Hz)", &draft_.f1, 0, 0, "%.4g");
    if (band) {
        ImGui::SetNextItemWidth(fh * 8);
        ImGui::InputDouble("upper corner (Hz)", &draft_.f2, 0, 0, "%.4g");
    }
    ImGui::Checkbox("zero phase (forward and backward, no time shift)", &draft_.zero_phase);
    // Check against the slowest and fastest sample rates in the data.
    for (const Trace& t : traces_) {
        std::string why = check_filter(draft_, t.rate);
        if (!why.empty()) {
            ImGui::TextColored(ImVec4(1, 0.6f, 0.3f, 1), "%s (%.4g samples/s): %s", t.id.c_str(), t.rate, why.c_str());
            break;
        }
    }
    if (!filters_.empty() && ImGui::Button("Save changes")) {
        filters_[edit_filter_] = draft_;
        save_filters(filters_);
    }
    if (!filters_.empty()) ImGui::SameLine();
    if (ImGui::Button("Save as new filter")) {
        filters_.push_back(draft_);
        edit_filter_ = (int)filters_.size() - 1;
        draft_src_ = edit_filter_;
        save_filters(filters_);
    }
    if (!filters_.empty()) {
        ImGui::SameLine();
        if (ImGui::Button("Delete")) {
            int d = edit_filter_;
            filters_.erase(filters_.begin() + d);
            auto fix = [&](int& f) {
                if (f == d) f = -1;
                else if (f > d) f--;
            };
            fix(proc_.filter);
            for (TraceOpts& to : topts_) fix(to.o.filter);
            draft_src_ = -2;
            save_filters(filters_);
        }
    }
    ImGui::EndGroup();
    ImGui::Separator();
    if (ImGui::Button("Close")) {
        save_filters(filters_);
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}

// ---------------------------------------------------------------------------
// Calibration window

void App::calib_editor() {
    if (!ImGui::BeginPopupModal("Calibration", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;
    const float fh = ImGui::GetFontSize();
    calib_sel_.resize(traces_.size(), 0);
    ImGui::PushTextWrapPos(fh * 42);
    ImGui::TextWrapped(
        "Bitweight is the number of physical units per count (for example 2.5e-6 for Pa per count). Once set, the "
        "Main, Magnify, Spectra and Spectrogram tabs show physical units. Values are remembered by channel code.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();
    ImGui::SetNextItemWidth(fh * 8);
    ImGui::InputText("##match", calib_match_, sizeof calib_match_);
    ImGui::SameLine();
    if (ImGui::Button("Tick channels containing this text"))
        for (size_t i = 0; i < traces_.size(); i++)
            calib_sel_[i] = *calib_match_ && traces_[i].id.find(calib_match_) != std::string::npos;
    ImGui::SameLine();
    if (ImGui::Button("Tick all")) std::fill(calib_sel_.begin(), calib_sel_.end(), 1);
    ImGui::SameLine();
    if (ImGui::Button("Tick none")) std::fill(calib_sel_.begin(), calib_sel_.end(), 0);

    ImGui::SetNextItemWidth(fh * 9);
    ImGui::InputDouble("bitweight##set", &calib_set_scale_, 0, 0, "%.6g");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(fh * 6);
    ImGui::InputText("units##set", calib_set_units_, sizeof calib_set_units_);
    ImGui::SameLine();
    if (ImGui::Button("Apply to ticked channels") && calib_set_scale_ != 0) {
        for (size_t i = 0; i < traces_.size(); i++)
            if (calib_sel_[i]) {
                Calib c;
                c.scale = calib_set_scale_;
                c.units = *calib_set_units_ ? calib_set_units_ : "counts";
                calib_[traces_[i].id] = c;
            }
    }
    ImGui::Spacing();
    if (ImGui::BeginTable("calib", 4, ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg | ImGuiTableFlags_Borders,
                          ImVec2(fh * 42, fh * 14))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, fh * 2);
        ImGui::TableSetupColumn("channel", ImGuiTableColumnFlags_WidthFixed, fh * 12);
        ImGui::TableSetupColumn("bitweight", ImGuiTableColumnFlags_WidthFixed, fh * 14);
        ImGui::TableSetupColumn("units", ImGuiTableColumnFlags_WidthFixed, fh * 8);
        ImGui::TableHeadersRow();
        for (int k = 0; k < (int)order_.size(); k++) {
            int ti = order_[k];
            ImGui::PushID(ti);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            bool tick = calib_sel_[ti] != 0;
            if (ImGui::Checkbox("##t", &tick)) calib_sel_[ti] = tick;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(traces_[ti].id.c_str());
            Calib c = calib_for(ti);
            bool commit = false;
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-1);
            ImGui::InputDouble("##bw", &c.scale, 0, 0, "%.6g");
            commit |= ImGui::IsItemDeactivatedAfterEdit();
            ImGui::TableNextColumn();
            char u[32];
            std::snprintf(u, sizeof u, "%s", c.units.c_str());
            ImGui::SetNextItemWidth(-1);
            if (ImGui::InputText("##u", u, sizeof u)) c.units = u;
            commit |= ImGui::IsItemDeactivatedAfterEdit();
            if (commit && c.scale != 0) {
                if (c.units.empty()) c.units = "counts";
                calib_[traces_[ti].id] = c;
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    if (ImGui::Button("Set all back to counts")) calib_.clear();
    ImGui::SameLine();
    if (ImGui::Button("Close")) {
        save_calib(calib_);
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndPopup();
}

// ---------------------------------------------------------------------------
// Magnify

void App::mag_click(double x) {
    if (pick_state_ == 0) {
        pick_first_ = x;
        pick_state_ = 1;
        return;
    }
    pick_state_ = 0;
    double a = std::min(pick_first_, x), b = std::max(pick_first_, x);
    if (b - a < 1e-6) return;
    view_[0] = a;  // the processed window follows by itself (see grow_window_if_needed)
    view_[1] = b;
}

void App::magnify_toolbar(MagState& ms) {
    ImGui::PushID(&ms);
    ImGui::Checkbox("Overlay", &ms.overlay);
    ImGui::SameLine(0, 20);
    combo_filter("Filter", proc_.filter);
    ImGui::SameLine();
    if (ImGui::Button("Edit filters...")) {
        if (proc_.filter >= 0) edit_filter_ = proc_.filter;
        ImGui::OpenPopup("Filters");
    }
    filter_editor();
    ImGui::SameLine(0, 20);
    combo_detrend("Detrend", proc_.detrend);
    ImGui::SameLine();
    ImGui::Checkbox("Reverse polarity", &proc_.reverse);

    ImGui::TextUnformatted("Scale:");
    ImGui::SameLine();
    if (ImGui::RadioButton("each trace", !ms.global_scale)) ms.global_scale = false;
    ImGui::SameLine();
    if (ImGui::RadioButton("all traces together", ms.global_scale)) ms.global_scale = true;
    ImGui::SameLine();
    ImGui::Checkbox("fit to the part on screen", &ms.fit_visible);
    ImGui::SameLine(0, 20);
    if (pick_state_ == 1)
        ImGui::Text("New t1 = %s. Click again to set t2 and zoom.", fmt_disp(pick_first_).c_str());
    else
        ImGui::Text("Showing %s to %s (%.3f s). Click twice on a plot to zoom to a new window.",
                    fmt_disp(view_[0]).c_str(), fmt_disp(view_[1]).c_str(), view_[1] - view_[0]);

    // Moving statistics, drawn under each trace (not available when traces are overlaid).
    ImGui::BeginDisabled(ms.overlay);
    ImGui::TextUnformatted("Stats:");
    ImGui::SameLine();
    ImGui::Checkbox("moving RMS", &ms.st_rms);
    ImGui::SameLine();
    ImGui::Checkbox("skewness", &ms.st_skew);
    ImGui::SameLine();
    ImGui::Checkbox("kurtosis", &ms.st_kurt);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5);
    ImGui::InputDouble("window (s)", &ms.st_win, 0, 0, "%.3g");
    if (!(ms.st_win > 0)) ms.st_win = 5.0;
    ImGui::EndDisabled();
    if (ms.overlay && (ms.st_rms || ms.st_skew || ms.st_kurt)) {
        ImGui::SameLine();
        ImGui::TextDisabled("(not shown in overlay mode)");
    }
    ImGui::PopID();
}

void App::draw_magnify(const char* id, float height, MagState& ms) {
    ImGui::BeginChild(id, ImVec2(0, height), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    if (!have_win_) {
        need_window_message();
        ImGui::EndChild();
        return;
    }
    magnify_toolbar(ms);
    std::vector<int> rows = page_traces(true);
    if (rows.empty()) {
        ImGui::TextUnformatted("All traces on this page are turned off (see the Main tab).");
        ImGui::EndChild();
        return;
    }
    ensure_processed(rows, true);

    const float label_w = ImGui::GetFontSize() * 11.0f;
    const float avail_h = ImGui::GetContentRegionAvail().y;
    const int nbins = std::max(100, std::min(4000, (int)(ImGui::GetContentRegionAvail().x - label_w - 20)));
    const int N = ms.overlay ? 1 : (int)rows.size();
    const float row_h = std::max(30.0f, std::floor((avail_h - axis_label_height()) / N));

    // Points to draw, and the vertical range, for each row.
    std::vector<std::vector<Series>> data(rows.size());
    std::vector<std::string> notes(rows.size());
    std::vector<double> ylo(rows.size(), 1e300), yhi(rows.size(), -1e300);
    for (size_t r = 0; r < rows.size(); r++) {
        auto P = proc_for(rows[r], true);
        if (!P) continue;
        double o = off(rows[r]);
        for (const Segment& seg : P->segs) {
            Series ser;
            seg.envelope(o + view_[0], o + view_[1], nbins, ser);
            if (!ser.x.empty()) {
                for (size_t i = 0; i < ser.x.size(); i++) {
                    if (ms.fit_visible) {
                        ylo[r] = std::min(ylo[r], ser.y[i]);
                        yhi[r] = std::max(yhi[r], ser.y[i]);
                    }
                    ser.x[i] -= o;
                }
                data[r].push_back(std::move(ser));
            }
            if (!ms.fit_visible && seg.size() > 0) {
                double mn, mx;
                seg.minmax(0, seg.size(), mn, mx);
                ylo[r] = std::min(ylo[r], mn);
                yhi[r] = std::max(yhi[r], mx);
            }
        }
        double rms = rms_for(rows[r], true);
        const std::string& u = calib_for(rows[r]).units;
        char buf[128];
        if (std::isnan(rms))
            std::snprintf(buf, sizeof buf, "no data in view");
        else
            std::snprintf(buf, sizeof buf, "RMS %.4g %s", rms, u == "counts" ? "" : u.c_str());
        notes[r] = buf;
        if (!P->warning.empty()) notes[r] += "\n" + P->warning;
    }
    if (ms.global_scale && !ms.overlay) {
        // Traces share a scale only when they have the same units.
        std::map<std::string, std::pair<double, double>> common;
        for (size_t r = 0; r < rows.size(); r++) {
            auto it = common.emplace(calib_for(rows[r]).units, std::make_pair(1e300, -1e300)).first;
            it->second.first = std::min(it->second.first, ylo[r]);
            it->second.second = std::max(it->second.second, yhi[r]);
        }
        for (size_t r = 0; r < rows.size(); r++) {
            const auto& c = common[calib_for(rows[r]).units];
            ylo[r] = c.first;
            yhi[r] = c.second;
        }
    }
    for (size_t r = 0; r < rows.size(); r++) {
        if (ylo[r] > yhi[r]) {
            ylo[r] = -1;
            yhi[r] = 1;
        }
        double span = yhi[r] - ylo[r];
        double big = std::max(std::fabs(ylo[r]), std::fabs(yhi[r]));
        if (span <= 1e-12 * big || big == 0) {  // flat trace: give it some room
            double m = big > 0 ? big : 1.0;
            ylo[r] -= m;
            yhi[r] += m;
        } else {
            ylo[r] -= 0.05 * span;
            yhi[r] += 0.05 * span;
        }
    }

    double dlo, dhi;
    data_range(dlo, dhi);
    const double dspan = dhi - dlo;

    bool clicked = false;
    double click_x = 0;
    auto mark_and_click = [&]() {
        if (pick_state_ == 1) {
            ImPlotSpec mark;
            mark.LineColor = mark_color();
            mark.LineWeight = 1.5f;
            mark.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##pick", &pick_first_, 1, mark);
        }
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            ImVec2 d = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
            if (d.x == 0 && d.y == 0) {
                clicked = true;
                click_x = ImPlot::GetPlotMousePos().x;
            }
        }
    };
    const ImPlotFlags pflags =
        ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
    auto setup_x = [&](bool last_row) {
        ImPlot::SetupAxis(ImAxis_X1, nullptr,
                          ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight |
                              (last_row ? 0 : ImPlotAxisFlags_NoTickLabels));
        ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
        ImPlot::SetupAxisLinks(ImAxis_X1, &view_[0], &view_[1]);
        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, dlo - dspan, dhi + dspan);  // do not get lost
    };

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
    if (ms.overlay) {
        ImGui::BeginChild("lbls", ImVec2(label_w, row_h));
        for (size_t r = 0; r < rows.size(); r++) {
            ImGui::PushID(rows[r]);
            trace_label(rows[r], false, true, notes[r].c_str());
            ImGui::Spacing();
            ImGui::PopID();
        }
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        const std::string units = units_of(rows);
        if (ImPlot::BeginPlot("##overlay", ImVec2(-1, row_h + axis_label_height()), pflags)) {
            setup_x(true);
            ImPlot::SetupAxis(ImAxis_Y1, units.c_str(), ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_NoHighlight);
            for (size_t r = 0; r < rows.size(); r++) {
                ImPlotSpec spec;
                spec.LineColor = trace_color(rows[r]);
                spec.Flags = ImPlotItemFlags_NoLegend;
                for (const Series& s : data[r])
                    ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);
            }
            mark_and_click();
            ImPlot::EndPlot();
        }
    } else {
        // Each trace gets its plot, then a small plot below it for every statistic switched on.
        const int nstat = (int)ms.st_rms + (int)ms.st_skew + (int)ms.st_kurt;
        std::vector<std::shared_ptr<const StatsResult>> stats;
        if (nstat > 0) stats = stats_for(rows, ms, nbins, id);
        const float stat_h =
            nstat > 0 ? std::max(18.0f, std::min(0.25f * row_h, 0.6f * row_h / (float)nstat)) : 0.0f;
        const float trace_h = std::max(25.0f, row_h - nstat * stat_h);

        // Draws one statistic's plot (kind 0 = RMS, 1 = skewness, 2 = kurtosis).
        auto stat_plot = [&](int ti, size_t r, int kind, bool last) {
            static const char* names[3] = {"RMS", "skewness", "kurtosis"};
            const std::vector<Series>* ser = nullptr;
            if (stats[r]) ser = kind == 0 ? &stats[r]->rms : kind == 1 ? &stats[r]->skew : &stats[r]->kurt;
            double lo = 1e300, hi = -1e300;
            if (ser)
                for (const Series& s : *ser)
                    for (size_t i = 0; i < s.y.size(); i++) {
                        if (s.x[i] < off(ti) + view_[0] || s.x[i] > off(ti) + view_[1]) continue;
                        lo = std::min(lo, s.y[i]);
                        hi = std::max(hi, s.y[i]);
                    }
            const double ref = kind == 1 ? 0.0 : kind == 2 ? 3.0 : 0.0;  // value for a normal distribution
            if (kind > 0 && lo <= hi) {  // keep the reference value in view
                lo = std::min(lo, ref);
                hi = std::max(hi, ref);
            }
            if (lo > hi) {
                lo = 0;
                hi = 1;
            }
            double span = hi - lo, big = std::max(std::fabs(lo), std::fabs(hi));
            if (span <= 1e-12 * big || big == 0) {
                double m = big > 0 ? big : 1.0;
                lo -= m;
                hi += m;
            } else {
                lo -= 0.08 * span;
                hi += 0.08 * span;
            }
            ImGui::PushID(kind + 100);
            char lbl[64];
            std::snprintf(lbl, sizeof lbl, "%s (%.3g s)", names[kind], ms.st_win);
            ImGui::BeginChild("slbl", ImVec2(label_w, stat_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
            ImGui::TextColored(trace_color(ti), "%s", lbl);
            ImGui::EndChild();
            ImGui::SameLine(0, 0);
            if (ImPlot::BeginPlot("##s", ImVec2(-1, last ? stat_h + axis_label_height() : stat_h), pflags)) {
                setup_x(last);
                const std::string yl = kind == 0 ? calib_for(ti).units : std::string();
                ImPlot::SetupAxis(ImAxis_Y1, yl.empty() ? nullptr : yl.c_str(),
                                  ImPlotAxisFlags_NoHighlight | ImPlotAxisFlags_Lock | (yl.empty() ? ImPlotAxisFlags_NoLabel : 0));
                ImPlot::SetupAxisLimits(ImAxis_Y1, lo, hi, ImPlotCond_Always);
                if (kind > 0) {  // faint reference line
                    ImPlotSpec rs;
                    rs.LineColor = ImVec4(0.6f, 0.6f, 0.6f, 0.5f);
                    rs.Flags = ImPlotItemFlags_NoLegend | ImPlotInfLinesFlags_Horizontal;
                    double rv = ref;
                    ImPlot::PlotInfLines("##ref", &rv, 1, rs);
                }
                if (ser) {
                    ImPlotSpec spec;
                    spec.LineColor = trace_color(ti);
                    spec.Flags = ImPlotItemFlags_NoLegend;
                    for (const Series& s : *ser) {
                        std::vector<double> xs(s.x);
                        for (double& x : xs) x -= off(ti);
                        ImPlot::PlotLine("##l", xs.data(), s.y.data(), (int)xs.size(), spec);
                    }
                }
                mark_and_click();
                ImPlot::EndPlot();
            }
            ImGui::PopID();
        };

        bool aligned = ImPlot::BeginAlignedPlots("##align_mag");
        for (size_t r = 0; r < rows.size(); r++) {
            int ti = rows[r];
            ImGui::PushID(ti);
            ImGui::BeginChild("lbl", ImVec2(label_w, trace_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
            trace_label(ti, false, true, notes[r].c_str());
            ImGui::EndChild();
            ImGui::SameLine(0, 0);
            const bool last = (r + 1 == rows.size()) && nstat == 0;
            const std::string units = calib_for(ti).units;
            if (ImPlot::BeginPlot("##p", ImVec2(-1, last ? trace_h + axis_label_height() : trace_h), pflags)) {
                setup_x(last);
                ImPlot::SetupAxis(ImAxis_Y1, units.c_str(), ImPlotAxisFlags_NoHighlight | ImPlotAxisFlags_Lock);
                ImPlot::SetupAxisLimits(ImAxis_Y1, ylo[r], yhi[r], ImPlotCond_Always);
                ImPlotSpec spec;
                spec.LineColor = trace_color(ti);
                spec.Flags = ImPlotItemFlags_NoLegend;
                for (const Series& s : data[r])
                    ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);
                mark_and_click();
                ImPlot::EndPlot();
            }
            const bool final_row = (r + 1 == rows.size());
            int shown = 0;
            for (int kind = 0; kind < 3; kind++) {
                bool on = kind == 0 ? ms.st_rms : kind == 1 ? ms.st_skew : ms.st_kurt;
                if (!on) continue;
                shown++;
                stat_plot(ti, r, kind, final_row && shown == nstat);
            }
            ImGui::PopID();
        }
        if (aligned) ImPlot::EndAlignedPlots();
    }
    ImGui::PopStyleVar();
    ImGui::EndChild();
    if (clicked) mag_click(click_x);
}

// ---------------------------------------------------------------------------
// Spectra

void App::spectra_toolbar(SpecState& ss) {
    ImGui::PushID(&ss);
    ImGui::Checkbox("Overlay", &ss.overlay);
    ImGui::SameLine();
    ImGui::Checkbox("Normalize", &ss.normalize);
    ImGui::SameLine();
    ImGui::Checkbox("Log frequency", &ss.log_x);
    ImGui::SameLine();
    ImGui::Checkbox("Log power", &ss.log_y);
    ImGui::SameLine(0, 20);
    ImGui::Checkbox("Same filter as Magnify tab", &ss.same_filter);
    if (!ss.same_filter) {
        ImGui::SameLine();
        ImGui::Checkbox("Filtered", &ss.filtered);
    }
    ImGui::Checkbox("Welch average", &ss.welch);
    if (ss.welch) {
        ImGui::SameLine();
        combo_pow2("window (samples)", ss.nwin, 64, 65536);
    }
    ImGui::SameLine(0, 20);
    ImGui::TextUnformatted("Axis limits (0 = automatic):");
    const float w = ImGui::GetFontSize() * 5;
    const char* names[4] = {"f min", "f max", "power min", "power max"};
    for (int i = 0; i < 4; i++) {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(w);
        ImGui::InputDouble(names[i], &ss.lim[i], 0, 0, "%.3g");
    }
    ImGui::PopID();
}

void App::draw_spectra(const char* id, float height, SpecState& ss) {
    ImGui::BeginChild(id, ImVec2(0, height), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
    if (!have_win_) {
        need_window_message();
        ImGui::EndChild();
        return;
    }
    spectra_toolbar(ss);
    std::vector<int> rows = page_traces(true);
    if (rows.empty()) {
        ImGui::TextUnformatted("All traces on this page are turned off (see the Main tab).");
        ImGui::EndChild();
        return;
    }
    const bool wf = ss.same_filter ? true : ss.filtered;
    ensure_processed(rows, wf);
    auto sp = spectra_for(rows, wf, ss, id);

    // Shared frequency axis: reset when the choice of spectra or mode changes.
    {
        std::string sig;
        for (int ti : rows) sig += traces_[ti].id + ",";
        sig += std::to_string((int)wf) + std::to_string((int)ss.welch) + std::to_string(ss.nwin) +
               std::to_string((int)ss.log_x);
        double fmin = 1e300, fmax = -1e300;
        for (const auto& r : sp)
            if (r && !r->pts.x.empty()) {
                fmin = std::min(fmin, r->pts.x.front());
                fmax = std::max(fmax, r->pts.x.back());
            }
        double want[2];
        bool have = fmin <= fmax;
        if (ss.lim[0] > 0 && ss.lim[1] > ss.lim[0]) {
            want[0] = ss.lim[0];
            want[1] = ss.lim[1];
            have = true;
        } else {
            want[0] = fmin;
            want[1] = fmax;
        }
        if (have && (want[0] != ss.last_desired[0] || want[1] != ss.last_desired[1] || sig != ss.last_sig)) {
            ss.xv[0] = ss.last_desired[0] = want[0];
            ss.xv[1] = ss.last_desired[1] = want[1];
            ss.last_sig = sig;
        }
    }

    const float label_w = ImGui::GetFontSize() * 11.0f;
    const float avail_h = ImGui::GetContentRegionAvail().y;
    const int N = ss.overlay ? 1 : (int)rows.size();
    const float row_h = std::max(30.0f, std::floor((avail_h - axis_label_height()) / N));
    const bool user_y = ss.lim[2] > 0 && ss.lim[3] > ss.lim[2];
    const ImPlotFlags pflags = ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

    auto setup_axes = [&](bool last_row, const std::string& ylabel) {
        ImPlot::SetupAxis(ImAxis_X1, nullptr,
                          ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit | ImPlotAxisFlags_NoHighlight |
                              (last_row ? 0 : ImPlotAxisFlags_NoTickLabels));
        ImPlot::SetupAxis(ImAxis_Y1, ylabel.c_str(),
                          ImPlotAxisFlags_NoHighlight | (user_y ? 0 : ImPlotAxisFlags_AutoFit));
        ImPlot::SetupAxisScale(ImAxis_X1, ss.log_x ? ImPlotScale_Log10 : ImPlotScale_Linear);
        ImPlot::SetupAxisScale(ImAxis_Y1, ss.log_y ? ImPlotScale_Log10 : ImPlotScale_Linear);
        ImPlot::SetupAxisLinks(ImAxis_X1, &ss.xv[0], &ss.xv[1]);
        if (user_y) ImPlot::SetupAxisLimits(ImAxis_Y1, ss.lim[2], ss.lim[3], ImPlotCond_Always);
    };
    // Power axis label: units squared per Hz (or "relative power" when normalized).
    auto spec_label = [&](const std::vector<int>& rs) {
        if (ss.normalize) return std::string("relative power");
        std::string u = units_of(rs);
        return u == "mixed units" ? u : u + "\xc2\xb2/Hz";
    };
    auto draw_one = [&](int ti, const SpecResult& r) {
        ImPlotSpec spec;
        spec.LineColor = trace_color(ti);
        spec.Flags = ImPlotItemFlags_NoLegend;
        if (ss.normalize && r.pmax > 0) {
            std::vector<double> y(r.pts.y.size());
            for (size_t i = 0; i < y.size(); i++) y[i] = r.pts.y[i] / r.pmax;
            ImPlot::PlotLine("##l", r.pts.x.data(), y.data(), (int)y.size(), spec);
        } else {
            ImPlot::PlotLine("##l", r.pts.x.data(), r.pts.y.data(), (int)r.pts.y.size(), spec);
        }
    };

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
    if (ss.overlay) {
        ImGui::BeginChild("lbls", ImVec2(label_w, row_h));
        for (int ti : rows) {
            ImGui::PushID(ti);
            trace_label(ti, false, false, nullptr);
            ImGui::PopID();
        }
        ImGui::EndChild();
        ImGui::SameLine(0, 0);
        if (ImPlot::BeginPlot("##overlay", ImVec2(-1, row_h + axis_label_height()), pflags)) {
            setup_axes(true, spec_label(rows));
            for (size_t r = 0; r < rows.size(); r++)
                if (sp[r] && !sp[r]->pts.x.empty()) draw_one(rows[r], *sp[r]);
            ImPlot::EndPlot();
        }
    } else {
        bool aligned = ImPlot::BeginAlignedPlots("##align_spec");
        for (size_t r = 0; r < rows.size(); r++) {
            int ti = rows[r];
            ImGui::PushID(ti);
            ImGui::BeginChild("lbl", ImVec2(label_w, row_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
            trace_label(ti, false, false, (sp[r] && sp[r]->pts.x.empty()) ? "no signal in view" : nullptr);
            ImGui::EndChild();
            ImGui::SameLine(0, 0);
            if (ImPlot::BeginPlot("##p", ImVec2(-1, (r + 1 == rows.size()) ? row_h + axis_label_height() : row_h), pflags)) {
                setup_axes(r + 1 == rows.size(), spec_label({ti}));
                if (sp[r] && !sp[r]->pts.x.empty()) draw_one(ti, *sp[r]);
                ImPlot::EndPlot();
            }
            ImGui::PopID();
        }
        if (aligned) ImPlot::EndAlignedPlots();
    }
    ImGui::PopStyleVar();
    ImGui::EndChild();
}

// ---------------------------------------------------------------------------

void App::draw_split() {
    if (!have_win_) {
        need_window_message();
        return;
    }
    float h = ImGui::GetContentRegionAvail().y;
    draw_spectra("##split_spec", std::floor(h * 0.45f), split_spec_);
    draw_magnify("##split_mag", ImGui::GetContentRegionAvail().y, split_mag_);
}

}  // namespace msv
