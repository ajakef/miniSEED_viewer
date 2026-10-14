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
    win_stack_.push_back({win_[0], win_[1]});
    apply_window(a, b);
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

    ImGui::BeginDisabled(win_stack_.empty());
    if (ImGui::Button("Back")) {
        auto w = win_stack_.back();
        win_stack_.pop_back();
        apply_window(w[0], w[1]);
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Reset zoom")) {
        view_[0] = win_[0];
        view_[1] = win_[1];
    }
    ImGui::SameLine(0, 20);
    if (pick_state_ == 1)
        ImGui::Text("New t1 = %s. Click again to set t2 and zoom.", fmt_disp(pick_first_).c_str());
    else
        ImGui::Text("Window %s to %s (%.3f s). Click twice on a plot to zoom to a new window.",
                    fmt_disp(win_[0]).c_str(), fmt_disp(win_[1]).c_str(), win_[1] - win_[0]);
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

    // Points to draw for each row.
    std::vector<std::vector<Series>> data(rows.size());
    std::vector<std::string> notes(rows.size());
    for (size_t r = 0; r < rows.size(); r++) {
        auto P = proc_for(rows[r], true);
        if (!P) continue;
        double o = off(rows[r]);
        for (const Segment& seg : P->segs) {
            Series ser;
            seg.envelope(o + view_[0], o + view_[1], nbins, ser);
            if (ser.x.empty()) continue;
            for (double& x : ser.x) x -= o;
            data[r].push_back(std::move(ser));
        }
        double rms = rms_for(rows[r], true);
        char buf[96];
        if (std::isnan(rms))
            std::snprintf(buf, sizeof buf, "no data in view");
        else
            std::snprintf(buf, sizeof buf, "RMS %.4g", rms);
        notes[r] = buf;
        if (!P->warning.empty()) notes[r] += "\n" + P->warning;
    }

    bool clicked = false;
    double click_x = 0;
    auto mark_and_click = [&](bool hovered_ok) {
        if (pick_state_ == 1) {
            ImPlotSpec mark;
            mark.LineColor = ImVec4(1.0f, 0.35f, 0.35f, 1);
            mark.LineWeight = 1.5f;
            mark.Flags = ImPlotItemFlags_NoLegend;
            ImPlot::PlotInfLines("##pick", &pick_first_, 1, mark);
        }
        if (hovered_ok && ImPlot::IsPlotHovered() && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            ImVec2 d = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left);
            if (d.x == 0 && d.y == 0) {
                clicked = true;
                click_x = ImPlot::GetPlotMousePos().x;
            }
        }
    };
    const ImPlotFlags pflags =
        ImPlotFlags_NoTitle | ImPlotFlags_NoLegend | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

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
        if (ImPlot::BeginPlot("##overlay", ImVec2(-1, row_h + axis_label_height()), pflags)) {
            ImPlot::SetupAxis(ImAxis_X1, nullptr, ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit);
            ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
            ImPlot::SetupAxisLinks(ImAxis_X1, &view_[0], &view_[1]);
            for (size_t r = 0; r < rows.size(); r++) {
                ImPlotSpec spec;
                spec.LineColor = ImPlot::GetColormapColor(rows[r]);
                spec.Flags = ImPlotItemFlags_NoLegend;
                for (const Series& s : data[r])
                    ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);
            }
            mark_and_click(true);
            ImPlot::EndPlot();
        }
    } else {
        bool aligned = ImPlot::BeginAlignedPlots("##align_mag");
        for (size_t r = 0; r < rows.size(); r++) {
            int ti = rows[r];
            ImGui::PushID(ti);
            ImGui::BeginChild("lbl", ImVec2(label_w, row_h), ImGuiChildFlags_None, ImGuiWindowFlags_NoScrollbar);
            trace_label(ti, false, true, notes[r].c_str());
            ImGui::EndChild();
            ImGui::SameLine(0, 0);
            const bool last = (r + 1 == rows.size());
            if (ImPlot::BeginPlot("##p", ImVec2(-1, (r + 1 == rows.size()) ? row_h + axis_label_height() : row_h), pflags)) {
                ImPlot::SetupAxis(ImAxis_X1, nullptr,
                                  ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit |
                                      (last ? 0 : ImPlotAxisFlags_NoTickLabels));
                ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_X1, cfg.relative_time ? ImPlotScale_Linear : ImPlotScale_Time);
                ImPlot::SetupAxisLinks(ImAxis_X1, &view_[0], &view_[1]);
                ImPlotSpec spec;
                spec.LineColor = ImPlot::GetColormapColor(ti);
                spec.Flags = ImPlotItemFlags_NoLegend;
                for (const Series& s : data[r])
                    ImPlot::PlotLine("##l", s.x.data(), s.y.data(), (int)s.x.size(), spec);
                mark_and_click(true);
                ImPlot::EndPlot();
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
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7);
        ImGui::InputInt("window (samples)", &ss.nwin, 0, 0);
        ss.nwin = std::max(16, ss.nwin);
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

    auto setup_axes = [&](bool last_row) {
        ImPlot::SetupAxis(ImAxis_X1, nullptr,
                          ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoInitialFit |
                              (last_row ? 0 : ImPlotAxisFlags_NoTickLabels));
        ImPlot::SetupAxis(ImAxis_Y1, nullptr, ImPlotAxisFlags_NoLabel | (user_y ? 0 : ImPlotAxisFlags_AutoFit));
        ImPlot::SetupAxisScale(ImAxis_X1, ss.log_x ? ImPlotScale_Log10 : ImPlotScale_Linear);
        ImPlot::SetupAxisScale(ImAxis_Y1, ss.log_y ? ImPlotScale_Log10 : ImPlotScale_Linear);
        ImPlot::SetupAxisLinks(ImAxis_X1, &ss.xv[0], &ss.xv[1]);
        if (user_y) ImPlot::SetupAxisLimits(ImAxis_Y1, ss.lim[2], ss.lim[3], ImPlotCond_Always);
    };
    auto draw_one = [&](int ti, const SpecResult& r) {
        ImPlotSpec spec;
        spec.LineColor = ImPlot::GetColormapColor(ti);
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
            setup_axes(true);
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
                setup_axes(r + 1 == rows.size());
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
