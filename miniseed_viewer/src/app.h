// The user interface: everything drawn inside the program's window.
#pragma once
#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "data.h"
#include "imgui.h"
#include "dsp.h"
#include "picks.h"
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

// A trace after calibration/detrending/filtering, limited to the processed window.
struct Processed {
    std::vector<Segment> segs;
    std::string warning;
};

// A computed spectrum, already thinned out for drawing.
struct SpecResult {
    Series pts;  // x = Hz, y = power
    double pmax = 0;
};

// A computed spectrogram: power in dB on a regular grid, plus the picture made from it.
struct SgData {
    int w = 0, h = 0;           // columns (time) and rows (frequency, lowest first)
    std::vector<float> db;      // h * w values, NaN where there is no data
    double t0 = 0, t1 = 1;      // display time of the left and right edge
    double f0 = 0, f1 = 1;      // frequency at the bottom and top edge
    double p_lo = 0, p_hi = 0;  // typical low and high power (for automatic colours)
    double p_max = 0;           // the highest value anywhere
    unsigned tex = 0;           // OpenGL texture holding the picture (0 = none yet)
    std::string tex_sig;        // colours the picture was made with
    SgData() = default;
    SgData(const SgData&) = delete;
    SgData& operator=(const SgData&) = delete;
    ~SgData();
};

// Moving statistics for one trace: one series of points per gap-free piece.
struct StatsResult {
    std::vector<Series> rms, skew, kurt;  // x = absolute time
};

struct MagState {
    bool overlay = false;
    bool pick_mode = false;      // clicks add/remove event picks instead of choosing a zoom window
    bool st_rms = false, st_skew = false, st_kurt = false;  // moving statistics to show (not in overlay mode)
    double st_win = 5.0;                                    // their window length, seconds
    bool global_scale = false;  // scale all traces together instead of each on its own
    bool fit_visible = true;    // scale to what is on screen, not the whole processed window
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

struct SgState {
    int nwin = 1024;
    int overlap_pct = 50;
    bool log_f = true;
    bool same_filter = true;
    bool filtered = true;
    bool normalize = false;      // each trace's peak = top of the color scale
    double db_rel[2] = {-60, 0}; // color range (dB below each trace's peak) when normalized and not automatic
    int cmap = 0;                // index into the list of colour maps
    bool auto_db = true;
    double db[2] = {-100, -20};  // colour range when not automatic
    double flim[2] = {0, 0};     // frequency limits, 0 = automatic
    // bookkeeping for the shared frequency axis
    double yv[2] = {0.1, 100};
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
    bool loading_pending_ = false;   // a page of traces is about to be read from disk
    std::vector<long> last_used_;    // per trace: frame number when it was last on screen
    long frame_no_ = 0;

    std::vector<FilterSpec> filters_;
    ProcOpts proc_;
    struct TraceOpts {
        bool custom = false;
        ProcOpts o;
    };
    std::vector<TraceOpts> topts_;
    int edit_filter_ = 0;
    FilterSpec draft_;    // filter being edited in the filter window
    int draft_src_ = -2;  // which saved filter draft_ was copied from

    std::set<std::string> flips_;  // channels whose polarity is flipped (the same flip whatever filter is used)

    // event picks
    std::vector<Pick> picks_;
    std::string picks_path_;
    char picks_path_buf_[512] = "";
    bool picks_autosave_ = true;
    char pick_label_[32] = "P";
    bool show_picks_ = false;
    std::string picks_msg_;

    // calibration (bitweight) by channel code
    std::map<std::string, Calib> calib_;
    std::vector<char> calib_sel_;  // per trace: ticked in the calibration window
    char calib_match_[64] = "";
    double calib_set_scale_ = 1.0;
    char calib_set_units_[32] = "Pa";

    // selections, in display time (clock seconds, or seconds since the trace started)
    double main_lim_[2] = {0, 1};
    double sel_[2] = {0, 0};
    int sel_state_ = 0;  // 0 nothing, 1 first point set, 2 both set
    bool have_win_ = false;
    double view_[2] = {0, 0};  // what the magnify axes currently show
    double win_[2] = {0, 0};   // the stretch of time that has been processed (always covers view_)
    double applied_sel_[2] = {0, 0};  // the Main selection that the magnify view was last set from
    bool sel_applied_ = false;
    int pick_state_ = 0;
    double pick_first_ = 0;
    double last_view_[2] = {0, 0};
    double view_changed_ = 0;

    MagState mag_, split_mag_;
    SpecState spec_, split_spec_;
    SgState sg_;

    std::map<std::string, std::shared_ptr<const Processed>> proc_cache_;
    std::map<std::string, std::shared_ptr<const SpecResult>> spec_cache_;
    std::map<std::string, std::shared_ptr<const SpecResult>> spec_stale_;
    std::map<std::string, std::shared_ptr<SgData>> sg_cache_;
    std::map<std::string, std::shared_ptr<SgData>> sg_stale_;
    std::map<std::string, double> rms_cache_;
    std::map<std::string, std::shared_ptr<const StatsResult>> stats_cache_;
    std::map<std::string, std::shared_ptr<const StatsResult>> stats_stale_;

    void rebuild_order();
    void reset_views();
    void clear_caches();
    double off(int ti) const;
    void data_range(double& lo, double& hi) const;
    std::string fmt_disp(double d) const;
    int n_pages() const;
    std::vector<int> page_traces(bool only_on) const;
    const ProcOpts& opts_for(int ti) const;
    const Calib& calib_for(int ti) const;
    std::string units_of(const std::vector<int>& rows) const;
    bool view_stable() const;
    bool effective_reverse(int ti) const;
    void picks_changed();
    void draw_picks_window();

    // ---- processing (app.cpp)
    std::string proc_key(int ti, bool with_filter) const;
    std::shared_ptr<const Processed> compute_processed(int ti, bool with_filter) const;
    void ensure_processed(const std::vector<int>& rows, bool with_filter);
    std::shared_ptr<const Processed> proc_for(int ti, bool with_filter) const;
    double rms_for(int ti, bool with_filter);
    std::vector<std::shared_ptr<const SpecResult>> spectra_for(const std::vector<int>& rows, bool with_filter,
                                                               const SpecState& ss, const char* panel);
    std::vector<std::shared_ptr<const StatsResult>> stats_for(const std::vector<int>& rows, const MagState& ms,
                                                              int nbins, const char* panel);
    std::vector<std::shared_ptr<SgData>> specgrams_for(const std::vector<int>& rows, bool with_filter,
                                                       const SgState& sg);

    // ---- drawing (app.cpp: top, main, config; app_views.cpp: magnify, spectra, split, pop-up windows;
    //      app_specgram.cpp: spectrogram)
    void draw_top_bar();
    void draw_main();
    void draw_config();
    void draw_magnify(const char* id, float height, MagState& ms);
    void draw_spectra(const char* id, float height, SpecState& ss);
    void draw_specgram(const char* id, float height, SgState& sg);
    void draw_split();
    void magnify_toolbar(MagState& ms);
    void spectra_toolbar(SpecState& ss);
    void specgram_toolbar(SgState& sg);
    void filter_editor();
    void calib_editor();
    void combo_filter(const char* label, int& idx);
    void combo_detrend(const char* label, int& d);
    void trace_label(int ti, bool checkbox, bool settings_button, const char* line2);
    void need_window_message();
    void main_click(double x);
    void mag_click(double x);
    void apply_window(double a, double b);
    void apply_main_selection();
    ImVec4 trace_color(int ti) const;
    ImVec4 mark_color() const;
    bool combo_pow2(const char* label, int& v, int lo, int hi);
    void grow_window_if_needed();
    bool ensure_page_loaded();
};

}  // namespace msv
