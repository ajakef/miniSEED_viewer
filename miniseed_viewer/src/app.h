// The user interface: everything drawn inside the program's window.
#pragma once
#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "data.h"
#include "dsp.h"
#include "settings.h"

namespace msv {

// Extra height the bottom plot of a stack needs for its axis labels.
float axis_label_height();

// How traces are cleaned up before plotting in the magnify and spectra tabs.
struct ProcOpts {
    int filter = -1;  // index into the saved filters, -1 = none
    int detrend = (int)Detrend::Mean;
    bool reverse = false;
};

// A trace after detrending/filtering, limited to the chosen time window.
struct Processed {
    std::vector<Segment> segs;
    std::string warning;
};

// A computed spectrum, already thinned out for drawing.
struct SpecResult {
    Series pts;  // x = Hz, y = power
    double pmax = 0;
};

struct MagState {
    bool overlay = false;
};

struct SpecState {
    bool overlay = false;
    bool normalize = false;
    bool same_filter = true;  // use the magnify tab's filter
    bool filtered = true;     // only when same_filter is off
    bool log_x = true, log_y = true;
    bool welch = false;
    int nwin = 1024;
    double lim[4] = {0, 0, 0, 0};  // x min, x max, y min, y max; 0 = automatic
    // bookkeeping for the shared frequency axis
    double xv[2] = {0.1, 100};
    double last_desired[2] = {0, 0};
    std::string last_sig;
};

class App {
public:
    explicit App(const std::vector<std::string>& files);
    void add_files(const std::vector<std::string>& files);
    void frame();      // draws one frame (call inside a full-window ImGui window)
    void save_all();   // writes settings to disk
    bool busy() const; // true while the display is still changing and needs quick redraws
    Config cfg;

private:
    // ---- data and bookkeeping (app.cpp)
    std::vector<Trace> traces_;
    std::vector<int> order_;  // display order: indices into traces_
    std::vector<std::string> messages_;
    int n_files_ = 0;
    int page_ = 0;
    int goto_tab_ = -1;
    bool show_config_ = false;

    std::vector<FilterSpec> filters_;
    ProcOpts proc_;
    struct TraceOpts {
        bool custom = false;
        ProcOpts o;
    };
    std::vector<TraceOpts> topts_;
    int edit_filter_ = 0;
    FilterSpec draft_;       // filter being edited in the filter window
    int draft_src_ = -2;     // which saved filter draft_ was copied from

    // selections, in display time (clock seconds, or seconds since the trace started)
    double main_lim_[2] = {0, 1};
    double sel_[2] = {0, 0};
    int sel_state_ = 0;  // 0 nothing, 1 first point set, 2 both set
    bool have_win_ = false;
    double win_[2] = {0, 0};   // the stretch of time being processed
    double view_[2] = {0, 0};  // what the magnify axes currently show
    std::vector<std::array<double, 2>> win_stack_;
    int pick_state_ = 0;
    double pick_first_ = 0;
    double last_view_[2] = {0, 0};
    double view_changed_ = 0;

    MagState mag_, split_mag_;
    SpecState spec_, split_spec_;

    std::map<std::string, std::shared_ptr<const Processed>> proc_cache_;
    std::map<std::string, std::shared_ptr<const SpecResult>> spec_cache_;
    std::map<std::string, std::shared_ptr<const SpecResult>> spec_stale_;
    std::map<std::string, double> rms_cache_;

    void rebuild_order();
    void reset_views();
    void clear_caches();
    double off(int ti) const;
    void data_range(double& lo, double& hi) const;
    std::string fmt_disp(double d) const;
    int n_pages() const;
    std::vector<int> page_traces(bool only_on) const;
    const ProcOpts& opts_for(int ti) const;
    bool view_stable() const;

    // ---- processing (app.cpp)
    std::string proc_key(int ti, bool with_filter) const;
    std::shared_ptr<const Processed> compute_processed(int ti, bool with_filter) const;
    void ensure_processed(const std::vector<int>& rows, bool with_filter);
    std::shared_ptr<const Processed> proc_for(int ti, bool with_filter) const;
    double rms_for(int ti, bool with_filter);
    std::vector<std::shared_ptr<const SpecResult>> spectra_for(const std::vector<int>& rows, bool with_filter,
                                                               const SpecState& ss, const char* panel);

    // ---- drawing (app.cpp: top, main, config; app_views.cpp: the rest)
    void draw_top_bar();
    void draw_main();
    void draw_config();
    void draw_magnify(const char* id, float height, MagState& ms);
    void draw_spectra(const char* id, float height, SpecState& ss);
    void draw_split();
    void magnify_toolbar(MagState& ms);
    void spectra_toolbar(SpecState& ss);
    void filter_editor();
    void combo_filter(const char* label, int& idx);
    void combo_detrend(const char* label, int& d);
    void trace_label(int ti, bool checkbox, bool settings_button, const char* line2);
    void need_window_message();
    void main_click(double x);
    void mag_click(double x);
    void apply_window(double a, double b);
};

}  // namespace msv
