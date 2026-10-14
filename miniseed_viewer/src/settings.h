// Options that are remembered between runs, and the list of saved filters.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "data.h"
#include "dsp.h"

namespace msv {

enum class SortKey { CommandLine = 0, FileName, StartTime, Network, Station, Location, Channel };
const int kNumSortKeys = 7;
const char* sort_key_name(SortKey k);

struct SortItem {
    SortKey key;
    bool use;        // part of the sort order
    bool descending;
};

struct Config {
    int n_per_screen = 12;
    bool main_global_norm = false;     // main tab: scale all traces together
    bool main_norm_visible = false;    // main tab: scale to the part on screen, not the whole trace
    bool relative_time = false;        // seconds since start of each trace, not clock time (UTC)
    int detrend = (int)Detrend::Mean;  // default for the magnify tab
    bool reverse_polarity = false;
    bool welch = false;                // spectra: Welch average instead of one spectrum
    int welch_window = 1024;
    bool spec_log_x = true, spec_log_y = true;
    bool spec_normalize = false;
    bool overlay_magnify = false, overlay_spectra = false;
    bool split_overlay_magnify = false, split_overlay_spectra = false;
    float ui_scale = 1.0f;
    std::vector<SortItem> sort;        // most important first
    Config();
};

// Folder where settings are kept (created if needed).
std::string config_dir();

void load_config(Config& c);
void save_config(const Config& c);

// Physical units per count for one channel, and the name of those units.
struct Calib {
    double scale = 1.0;
    std::string units = "counts";
};
// Calibrations are kept by channel code (NET.STA.LOC.CHA).
void load_calib(std::map<std::string, Calib>& calib);
void save_calib(const std::map<std::string, Calib>& calib);

void load_filters(std::vector<FilterSpec>& filters);
void save_filters(const std::vector<FilterSpec>& filters);

// Display order for the traces: indices into traces, sorted by the keys in
// order (most important first). Ties keep their original order.
std::vector<int> sort_order(const std::vector<Trace>& traces, const std::vector<SortItem>& order);

}  // namespace msv
